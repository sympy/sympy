"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.3 Tangent/4.3.7 (d trig)^m (a+b (c tan)^n)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n, p = symbols('a b c d e f m n p')

def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_1():
    f = (b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -b**2*sqrt(b*tan(e + f*x)**2)*log(cos(e + f*x))*cot(e + f*x)/f + b**2*sqrt(b*tan(e + f*x)**2)*tan(e + f*x)**3/(4*f) - b**2*sqrt(b*tan(e + f*x)**2)*tan(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_2():
    f = (b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = b*sqrt(b*tan(e + f*x)**2)*log(cos(e + f*x))*cot(e + f*x)/f + b*sqrt(b*tan(e + f*x)**2)*tan(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_3():
    f = sqrt(b*tan(e + f*x)**2)
    F = -sqrt(b*tan(e + f*x)**2)*log(cos(e + f*x))*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_4():
    f = 1/sqrt(b*tan(e + f*x)**2)
    F = log(sin(e + f*x))*tan(e + f*x)/(f*sqrt(b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_5():
    f = (b*tan(e + f*x)**2)**(sympy.S(-3)/2)
    F = -log(sin(e + f*x))*tan(e + f*x)/(b*f*sqrt(b*tan(e + f*x)**2)) - cot(e + f*x)/(2*b*f*sqrt(b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_6():
    f = (b*tan(e + f*x)**2)**(sympy.S(-5)/2)
    F = log(sin(e + f*x))*tan(e + f*x)/(b**2*f*sqrt(b*tan(e + f*x)**2)) - cot(e + f*x)**3/(4*b**2*f*sqrt(b*tan(e + f*x)**2)) + cot(e + f*x)/(2*b**2*f*sqrt(b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_7():
    f = (b*tan(e + f*x)**3)**(sympy.S(5)/2)
    F = -sqrt(2)*b**2*sqrt(b*tan(e + f*x)**3)*log(-sqrt(2)*sqrt(tan(e + f*x)) + tan(e + f*x) + 1)/(4*f*tan(e + f*x)**(sympy.S(3)/2)) + sqrt(2)*b**2*sqrt(b*tan(e + f*x)**3)*log(sqrt(2)*sqrt(tan(e + f*x)) + tan(e + f*x) + 1)/(4*f*tan(e + f*x)**(sympy.S(3)/2)) + 2*b**2*sqrt(b*tan(e + f*x)**3)*tan(e + f*x)**5/(13*f) - 2*b**2*sqrt(b*tan(e + f*x)**3)*tan(e + f*x)**3/(9*f) + 2*b**2*sqrt(b*tan(e + f*x)**3)*tan(e + f*x)/(5*f) - 2*b**2*sqrt(b*tan(e + f*x)**3)*cot(e + f*x)/f + sqrt(2)*b**2*sqrt(b*tan(e + f*x)**3)*atan(sqrt(2)*sqrt(tan(e + f*x)) - 1)/(2*f*tan(e + f*x)**(sympy.S(3)/2)) + sqrt(2)*b**2*sqrt(b*tan(e + f*x)**3)*atan(sqrt(2)*sqrt(tan(e + f*x)) + 1)/(2*f*tan(e + f*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_8():
    f = (b*tan(e + f*x)**3)**(sympy.S(3)/2)
    F = sqrt(2)*b*sqrt(b*tan(e + f*x)**3)*log(-sqrt(2)*sqrt(tan(e + f*x)) + tan(e + f*x) + 1)/(4*f*tan(e + f*x)**(sympy.S(3)/2)) - sqrt(2)*b*sqrt(b*tan(e + f*x)**3)*log(sqrt(2)*sqrt(tan(e + f*x)) + tan(e + f*x) + 1)/(4*f*tan(e + f*x)**(sympy.S(3)/2)) + 2*b*sqrt(b*tan(e + f*x)**3)*tan(e + f*x)**2/(7*f) - 2*b*sqrt(b*tan(e + f*x)**3)/(3*f) + sqrt(2)*b*sqrt(b*tan(e + f*x)**3)*atan(sqrt(2)*sqrt(tan(e + f*x)) - 1)/(2*f*tan(e + f*x)**(sympy.S(3)/2)) + sqrt(2)*b*sqrt(b*tan(e + f*x)**3)*atan(sqrt(2)*sqrt(tan(e + f*x)) + 1)/(2*f*tan(e + f*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_9():
    f = sqrt(b*tan(e + f*x)**3)
    F = sqrt(2)*sqrt(b*tan(e + f*x)**3)*log(-sqrt(2)*sqrt(tan(e + f*x)) + tan(e + f*x) + 1)/(4*f*tan(e + f*x)**(sympy.S(3)/2)) - sqrt(2)*sqrt(b*tan(e + f*x)**3)*log(sqrt(2)*sqrt(tan(e + f*x)) + tan(e + f*x) + 1)/(4*f*tan(e + f*x)**(sympy.S(3)/2)) + 2*sqrt(b*tan(e + f*x)**3)*cot(e + f*x)/f - sqrt(2)*sqrt(b*tan(e + f*x)**3)*atan(sqrt(2)*sqrt(tan(e + f*x)) - 1)/(2*f*tan(e + f*x)**(sympy.S(3)/2)) - sqrt(2)*sqrt(b*tan(e + f*x)**3)*atan(sqrt(2)*sqrt(tan(e + f*x)) + 1)/(2*f*tan(e + f*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_10():
    f = 1/sqrt(b*tan(e + f*x)**3)
    F = -sqrt(2)*log(-sqrt(2)*sqrt(tan(e + f*x)) + tan(e + f*x) + 1)*tan(e + f*x)**(sympy.S(3)/2)/(4*f*sqrt(b*tan(e + f*x)**3)) + sqrt(2)*log(sqrt(2)*sqrt(tan(e + f*x)) + tan(e + f*x) + 1)*tan(e + f*x)**(sympy.S(3)/2)/(4*f*sqrt(b*tan(e + f*x)**3)) - sqrt(2)*tan(e + f*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(e + f*x)) - 1)/(2*f*sqrt(b*tan(e + f*x)**3)) - sqrt(2)*tan(e + f*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(e + f*x)) + 1)/(2*f*sqrt(b*tan(e + f*x)**3)) - 2*tan(e + f*x)/(f*sqrt(b*tan(e + f*x)**3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_11():
    f = (b*tan(e + f*x)**3)**(sympy.S(-3)/2)
    F = -sqrt(2)*log(-sqrt(2)*sqrt(tan(e + f*x)) + tan(e + f*x) + 1)*tan(e + f*x)**(sympy.S(3)/2)/(4*b*f*sqrt(b*tan(e + f*x)**3)) + sqrt(2)*log(sqrt(2)*sqrt(tan(e + f*x)) + tan(e + f*x) + 1)*tan(e + f*x)**(sympy.S(3)/2)/(4*b*f*sqrt(b*tan(e + f*x)**3)) + sqrt(2)*tan(e + f*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(e + f*x)) - 1)/(2*b*f*sqrt(b*tan(e + f*x)**3)) + sqrt(2)*tan(e + f*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(e + f*x)) + 1)/(2*b*f*sqrt(b*tan(e + f*x)**3)) - 2*cot(e + f*x)**2/(7*b*f*sqrt(b*tan(e + f*x)**3)) + 2/(3*b*f*sqrt(b*tan(e + f*x)**3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_12():
    f = (b*tan(e + f*x)**3)**(sympy.S(-5)/2)
    F = sqrt(2)*log(-sqrt(2)*sqrt(tan(e + f*x)) + tan(e + f*x) + 1)*tan(e + f*x)**(sympy.S(3)/2)/(4*b**2*f*sqrt(b*tan(e + f*x)**3)) - sqrt(2)*log(sqrt(2)*sqrt(tan(e + f*x)) + tan(e + f*x) + 1)*tan(e + f*x)**(sympy.S(3)/2)/(4*b**2*f*sqrt(b*tan(e + f*x)**3)) + sqrt(2)*tan(e + f*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(e + f*x)) - 1)/(2*b**2*f*sqrt(b*tan(e + f*x)**3)) + sqrt(2)*tan(e + f*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(e + f*x)) + 1)/(2*b**2*f*sqrt(b*tan(e + f*x)**3)) + 2*tan(e + f*x)/(b**2*f*sqrt(b*tan(e + f*x)**3)) - 2*cot(e + f*x)**5/(13*b**2*f*sqrt(b*tan(e + f*x)**3)) + 2*cot(e + f*x)**3/(9*b**2*f*sqrt(b*tan(e + f*x)**3)) - 2*cot(e + f*x)/(5*b**2*f*sqrt(b*tan(e + f*x)**3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_13():
    f = (b*tan(e + f*x)**4)**(sympy.S(5)/2)
    F = -b**2*x*sqrt(b*tan(e + f*x)**4)*cot(e + f*x)**2 + b**2*sqrt(b*tan(e + f*x)**4)*tan(e + f*x)**7/(9*f) - b**2*sqrt(b*tan(e + f*x)**4)*tan(e + f*x)**5/(7*f) + b**2*sqrt(b*tan(e + f*x)**4)*tan(e + f*x)**3/(5*f) - b**2*sqrt(b*tan(e + f*x)**4)*tan(e + f*x)/(3*f) + b**2*sqrt(b*tan(e + f*x)**4)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_14():
    f = (b*tan(e + f*x)**4)**(sympy.S(3)/2)
    F = -b*x*sqrt(b*tan(e + f*x)**4)*cot(e + f*x)**2 + b*sqrt(b*tan(e + f*x)**4)*tan(e + f*x)**3/(5*f) - b*sqrt(b*tan(e + f*x)**4)*tan(e + f*x)/(3*f) + b*sqrt(b*tan(e + f*x)**4)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_15():
    f = sqrt(b*tan(e + f*x)**4)
    F = -x*sqrt(b*tan(e + f*x)**4)*cot(e + f*x)**2 + sqrt(b*tan(e + f*x)**4)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_16():
    f = 1/sqrt(b*tan(e + f*x)**4)
    F = -x*tan(e + f*x)**2/sqrt(b*tan(e + f*x)**4) - tan(e + f*x)/(f*sqrt(b*tan(e + f*x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_17():
    f = (b*tan(e + f*x)**4)**(sympy.S(-3)/2)
    F = -x*tan(e + f*x)**2/(b*sqrt(b*tan(e + f*x)**4)) - tan(e + f*x)/(b*f*sqrt(b*tan(e + f*x)**4)) - cot(e + f*x)**3/(5*b*f*sqrt(b*tan(e + f*x)**4)) + cot(e + f*x)/(3*b*f*sqrt(b*tan(e + f*x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_18():
    f = (b*tan(e + f*x)**4)**(sympy.S(-5)/2)
    F = -x*tan(e + f*x)**2/(b**2*sqrt(b*tan(e + f*x)**4)) - tan(e + f*x)/(b**2*f*sqrt(b*tan(e + f*x)**4)) - cot(e + f*x)**7/(9*b**2*f*sqrt(b*tan(e + f*x)**4)) + cot(e + f*x)**5/(7*b**2*f*sqrt(b*tan(e + f*x)**4)) - cot(e + f*x)**3/(5*b**2*f*sqrt(b*tan(e + f*x)**4)) + cot(e + f*x)/(3*b**2*f*sqrt(b*tan(e + f*x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_19():
    f = (b*tan(e + f*x)**n)**(sympy.S(5)/2)
    F = 2*b**2*sqrt(b*tan(e + f*x)**n)*tan(e + f*x)**(2*n + 1)*hyper((1, 5*n/4 + sympy.S.Half), (5*n/4 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(5*n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_20():
    f = (b*tan(e + f*x)**n)**(sympy.S(3)/2)
    F = 2*b*sqrt(b*tan(e + f*x)**n)*tan(e + f*x)**(n + 1)*hyper((1, 3*n/4 + sympy.S.Half), (3*n/4 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(3*n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_21():
    f = sqrt(b*tan(e + f*x)**n)
    F = 2*sqrt(b*tan(e + f*x)**n)*tan(e + f*x)*hyper((1, n/4 + sympy.S.Half), (n/4 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_22():
    f = 1/sqrt(b*tan(e + f*x)**n)
    F = 2*tan(e + f*x)*hyper((1, sympy.S.Half - n/4), (sympy.S(3)/2 - n/4,), -tan(e + f*x)**2)/(f*sqrt(b*tan(e + f*x)**n)*(2 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_23():
    f = (b*tan(e + f*x)**n)**(sympy.S(-3)/2)
    F = 2*tan(e + f*x)**(1 - n)*hyper((1, sympy.S.Half - 3*n/4), (sympy.S(3)/2 - 3*n/4,), -tan(e + f*x)**2)/(b*f*sqrt(b*tan(e + f*x)**n)*(2 - 3*n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_24():
    f = (b*tan(e + f*x)**n)**(sympy.S(-5)/2)
    F = 2*tan(e + f*x)**(1 - 2*n)*hyper((1, sympy.S.Half - 5*n/4), (sympy.S(3)/2 - 5*n/4,), -tan(e + f*x)**2)/(b**2*f*sqrt(b*tan(e + f*x)**n)*(2 - 5*n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_25():
    f = (b*tan(e + f*x)**n)**p
    F = (b*tan(e + f*x)**n)**p*tan(e + f*x)*hyper((1, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_26():
    f = (b*tan(e + f*x)**2)**p
    F = (b*tan(e + f*x)**2)**p*tan(e + f*x)*hyper((1, p + sympy.S.Half), (p + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_27():
    f = (b*tan(e + f*x)**3)**p
    F = (b*tan(e + f*x)**3)**p*tan(e + f*x)*hyper((1, 3*p/2 + sympy.S.Half), (3*p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(3*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_28():
    f = (b*tan(e + f*x)**4)**p
    F = (b*tan(e + f*x)**4)**p*tan(e + f*x)*hyper((1, 2*p + sympy.S.Half), (2*p + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(4*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_29():
    f = (b*tan(e + f*x)**n)**(1/n)
    F = -(b*tan(e + f*x)**n)**(1/n)*log(cos(e + f*x))*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_30():
    f = (a + b*tan(e + f*x)**2)*sin(e + f*x)**5
    F = b*sec(e + f*x)/f - (a - 3*b)*cos(e + f*x)/f - (a - b)*cos(e + f*x)**5/(5*f) + (2*a - 3*b)*cos(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_31():
    f = (a + b*tan(e + f*x)**2)*sin(e + f*x)**3
    F = b*sec(e + f*x)/f - (a - 2*b)*cos(e + f*x)/f + (a - b)*cos(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_32():
    f = (a + b*tan(e + f*x)**2)*sin(e + f*x)
    F = b*sec(e + f*x)/f - (a - b)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_33():
    f = (a + b*tan(e + f*x)**2)*csc(e + f*x)
    F = -a*atanh(cos(e + f*x))/f + b*sec(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_34():
    f = (a + b*tan(e + f*x)**2)*csc(e + f*x)**3
    F = -a*cot(e + f*x)*csc(e + f*x)/(2*f) + b*sec(e + f*x)/f - (a + 2*b)*atanh(cos(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_35():
    f = (a + b*tan(e + f*x)**2)*csc(e + f*x)**5
    F = -a*cot(e + f*x)**3*csc(e + f*x)/(4*f) + b*sec(e + f*x)/f + (-3*a - 12*b)*atanh(cos(e + f*x))/(8*f) - (5*a + 4*b)*cot(e + f*x)*csc(e + f*x)/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_36():
    f = (a + b*tan(e + f*x)**2)*sin(e + f*x)**6
    F = b*tan(e + f*x)/f + x*(5*a - 35*b)/16 - (a - b)*sin(e + f*x)*cos(e + f*x)**5/(6*f) - (11*a - 29*b)*sin(e + f*x)*cos(e + f*x)/(16*f) + (13*a - 19*b)*sin(e + f*x)*cos(e + f*x)**3/(24*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_37():
    f = (a + b*tan(e + f*x)**2)*sin(e + f*x)**4
    F = b*tan(e + f*x)/f + x*(3*a - 15*b)/8 + (a - b)*sin(e + f*x)*cos(e + f*x)**3/(4*f) - (5*a - 9*b)*sin(e + f*x)*cos(e + f*x)/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_38():
    f = (a + b*tan(e + f*x)**2)*sin(e + f*x)**2
    F = b*tan(e + f*x)/f + x*(a - 3*b)/2 - (a - b)*sin(e + f*x)*cos(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_39():
    f = a + b*tan(e + f*x)**2
    F = a*x - b*x + b*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_40():
    f = (a + b*tan(e + f*x)**2)*csc(e + f*x)**2
    F = -a*cot(e + f*x)/f + b*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_41():
    f = (a + b*tan(e + f*x)**2)*csc(e + f*x)**4
    F = -a*cot(e + f*x)**3/(3*f) + b*tan(e + f*x)/f - (a + b)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_42():
    f = (a + b*tan(e + f*x)**2)*csc(e + f*x)**6
    F = -a*cot(e + f*x)**5/(5*f) + b*tan(e + f*x)/f - (a + 2*b)*cot(e + f*x)/f - (2*a + b)*cot(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_43():
    f = (a + b*tan(e + f*x)**2)**2*sin(e + f*x)**5
    F = b**2*sec(e + f*x)**3/(3*f) + b*(2*a - 4*b)*sec(e + f*x)/f - (a - b)**2*cos(e + f*x)**5/(5*f) + (a - b)*(2*a - 4*b)*cos(e + f*x)**3/(3*f) - (a**2 - 6*a*b + 6*b**2)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_44():
    f = (a + b*tan(e + f*x)**2)**2*sin(e + f*x)**3
    F = b**2*sec(e + f*x)**3/(3*f) + b*(2*a - 3*b)*sec(e + f*x)/f - (a - 3*b)*(a - b)*cos(e + f*x)/f + (a - b)**2*cos(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_45():
    f = (a + b*tan(e + f*x)**2)**2*sin(e + f*x)
    F = b**2*sec(e + f*x)**3/(3*f) + b*(2*a - 2*b)*sec(e + f*x)/f - (a - b)**2*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_46():
    f = (a + b*tan(e + f*x)**2)**2*csc(e + f*x)
    F = -a**2*atanh(cos(e + f*x))/f + b**2*sec(e + f*x)**3/(3*f) + b*(2*a - b)*sec(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_47():
    f = (a + b*tan(e + f*x)**2)**2*csc(e + f*x)**3
    F = -a**2*csc(e + f*x)**2*sec(e + f*x)/(2*f) - a*(a + 4*b)*atanh(cos(e + f*x))/(2*f) + a*(a + 4*b)*sec(e + f*x)/(2*f) + b**2*sec(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_48():
    f = (a + b*tan(e + f*x)**2)**2*csc(e + f*x)**5
    F = -a**2*csc(e + f*x)**4*sec(e + f*x)/(4*f) - a*(a + 8*b)*cot(e + f*x)*csc(e + f*x)/(8*f) + b**2*sec(e + f*x)**3/(3*f) + (a**2 + 8*a*b + 4*b**2)*sec(e + f*x)/(4*f) - (3*a**2 + 24*a*b + 8*b**2)*atanh(cos(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_49():
    f = (a + b*tan(e + f*x)**2)**2*sin(e + f*x)**4
    F = b**2*tan(e + f*x)**3/(3*f) + x*(3*a**2 - 30*a*b + 35*b**2)/8 - (a - 9*b)*(a - b)*sin(e + f*x)*cos(e + f*x)/(8*f) + (a - b)**2*sin(e + f*x)**4*tan(e + f*x)/(4*f) - (a**2 - 10*a*b + 13*b**2)*tan(e + f*x)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_50():
    f = (a + b*tan(e + f*x)**2)**2*sin(e + f*x)**2
    F = b**2*tan(e + f*x)**3/(3*f) + x*(a - 5*b)*(a - b)/2 - (a - 5*b)*(a - b)*tan(e + f*x)/(2*f) + (a - b)**2*sin(e + f*x)**2*tan(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_51():
    f = (a + b*tan(e + f*x)**2)**2
    F = b**2*tan(e + f*x)**3/(3*f) + b*(2*a - b)*tan(e + f*x)/f + x*(a - b)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_52():
    f = (a + b*tan(e + f*x)**2)**2*csc(e + f*x)**2
    F = -a**2*cot(e + f*x)/f + 2*a*b*tan(e + f*x)/f + b**2*tan(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_53():
    f = (a + b*tan(e + f*x)**2)**2*csc(e + f*x)**4
    F = -a**2*cot(e + f*x)**3/(3*f) - a*(a + 2*b)*cot(e + f*x)/f + b**2*tan(e + f*x)**3/(3*f) + b*(2*a + b)*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_54():
    f = (a + b*tan(e + f*x)**2)**2*csc(e + f*x)**6
    F = -a**2*cot(e + f*x)**5/(5*f) - 2*a*(a + b)*cot(e + f*x)**3/(3*f) + b**2*tan(e + f*x)**3/(3*f) + 2*b*(a + b)*tan(e + f*x)/f - (a**2 + 4*a*b + b**2)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_55():
    f = sin(e + f*x)**5/(a + b*tan(e + f*x)**2)
    F = -a**2*sqrt(b)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(f*(a - b)**(sympy.S(7)/2)) - a**2*cos(e + f*x)/(f*(a - b)**3) - cos(e + f*x)**5/(f*(5*a - 5*b)) + (2*a - b)*cos(e + f*x)**3/(3*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_56():
    f = sin(e + f*x)**3/(a + b*tan(e + f*x)**2)
    F = -a*sqrt(b)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(f*(a - b)**(sympy.S(5)/2)) - a*cos(e + f*x)/(f*(a - b)**2) + cos(e + f*x)**3/(f*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_57():
    f = sin(e + f*x)/(a + b*tan(e + f*x)**2)
    F = -sqrt(b)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(f*(a - b)**(sympy.S(3)/2)) - cos(e + f*x)/(f*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_58():
    f = csc(e + f*x)/(a + b*tan(e + f*x)**2)
    F = -sqrt(b)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(a*f*sqrt(a - b)) - atanh(cos(e + f*x))/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_59():
    f = csc(e + f*x)**3/(a + b*tan(e + f*x)**2)
    F = -cot(e + f*x)*csc(e + f*x)/(2*a*f) - sqrt(b)*sqrt(a - b)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(a**2*f) - (a - 2*b)*atanh(cos(e + f*x))/(2*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_60():
    f = csc(e + f*x)**5/(a + b*tan(e + f*x)**2)
    F = -cot(e + f*x)**3*csc(e + f*x)/(4*a*f) - (5*a - 4*b)*cot(e + f*x)*csc(e + f*x)/(8*a**2*f) - sqrt(b)*(a - b)**(sympy.S(3)/2)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(a**3*f) - (3*a**2 - 12*a*b + 8*b**2)*atanh(cos(e + f*x))/(8*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_61():
    f = sin(e + f*x)**6/(a + b*tan(e + f*x)**2)
    F = -a**(sympy.S(5)/2)*sqrt(b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(f*(a - b)**4) + x*(5*a**3 + 15*a**2*b - 5*a*b**2 + b**3)/(16*(a - b)**4) + sin(e + f*x)**3*cos(e + f*x)**3/(f*(6*a - 6*b)) + (3*a - b)*sin(e + f*x)*cos(e + f*x)**3/(8*f*(a - b)**2) - (11*a**2 - 4*a*b + b**2)*sin(e + f*x)*cos(e + f*x)/(16*f*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_62():
    f = sin(e + f*x)**4/(a + b*tan(e + f*x)**2)
    F = -a**(sympy.S(3)/2)*sqrt(b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(f*(a - b)**3) + x*(3*a**2 + 6*a*b - b**2)/(8*(a - b)**3) + sin(e + f*x)*cos(e + f*x)**3/(f*(4*a - 4*b)) - (5*a - b)*sin(e + f*x)*cos(e + f*x)/(8*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_63():
    f = sin(e + f*x)**2/(a + b*tan(e + f*x)**2)
    F = -sqrt(a)*sqrt(b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(f*(a - b)**2) + x*(a + b)/(2*(a - b)**2) - sin(e + f*x)*cos(e + f*x)/(f*(2*a - 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_64():
    f = 1/(a + b*tan(e + f*x)**2)
    F = x/(a - b) - sqrt(b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(sqrt(a)*f*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_65():
    f = csc(e + f*x)**2/(a + b*tan(e + f*x)**2)
    F = -cot(e + f*x)/(a*f) - sqrt(b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_66():
    f = csc(e + f*x)**4/(a + b*tan(e + f*x)**2)
    F = -cot(e + f*x)**3/(3*a*f) - (a - b)*cot(e + f*x)/(a**2*f) - sqrt(b)*(a - b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_67():
    f = csc(e + f*x)**6/(a + b*tan(e + f*x)**2)
    F = -cot(e + f*x)**5/(5*a*f) - (2*a - b)*cot(e + f*x)**3/(3*a**2*f) - (a - b)**2*cot(e + f*x)/(a**3*f) - sqrt(b)*(a - b)**2*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_68():
    f = sin(e + f*x)**5/(a + b*tan(e + f*x)**2)**2
    F = -a*sqrt(b)*(3*a + 4*b)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(2*f*(a - b)**(sympy.S(9)/2)) - b*(5*a**2 + 2*b**2)*sec(e + f*x)/(10*f*(a - b)**4*(a + b*sec(e + f*x)**2 - b)) - cos(e + f*x)**5/(f*(5*a - 5*b)*(a + b*sec(e + f*x)**2 - b)) + (10*a - 3*b)*cos(e + f*x)**3/(15*f*(a - b)**3) - (5*a**2 + 10*a*b - b**2)*cos(e + f*x)/(5*f*(a - b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_69():
    f = sin(e + f*x)**3/(a + b*tan(e + f*x)**2)**2
    F = -a*b*sec(e + f*x)/(2*f*(a - b)**3*(a + b*sec(e + f*x)**2 - b)) - sqrt(b)*(3*a + 2*b)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(2*f*(a - b)**(sympy.S(7)/2)) + cos(e + f*x)**3/(3*f*(a - b)**2) - (a + b)*cos(e + f*x)/(f*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_70():
    f = sin(e + f*x)/(a + b*tan(e + f*x)**2)**2
    F = -3*sqrt(b)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(2*f*(a - b)**(sympy.S(5)/2)) + cos(e + f*x)/(f*(2*a - 2*b)*(a + b*sec(e + f*x)**2 - b)) - 3*cos(e + f*x)/(2*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_71():
    f = csc(e + f*x)/(a + b*tan(e + f*x)**2)**2
    F = -b*sec(e + f*x)/(2*a*f*(a - b)*(a + b*sec(e + f*x)**2 - b)) - sqrt(b)*(3*a - 2*b)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(2*a**2*f*(a - b)**(sympy.S(3)/2)) - atanh(cos(e + f*x))/(a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_72():
    f = csc(e + f*x)**3/(a + b*tan(e + f*x)**2)**2
    F = -cot(e + f*x)*csc(e + f*x)/(2*a*f*(a + b*sec(e + f*x)**2 - b)) - b*sec(e + f*x)/(a**2*f*(a + b*sec(e + f*x)**2 - b)) - sqrt(b)*(3*a - 4*b)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(2*a**3*f*sqrt(a - b)) - (a - 4*b)*atanh(cos(e + f*x))/(2*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_73():
    f = csc(e + f*x)**5/(a + b*tan(e + f*x)**2)**2
    F = -cot(e + f*x)**3*csc(e + f*x)/(4*a*f*(a + b*sec(e + f*x)**2 - b)) - (5*a - 6*b)*cot(e + f*x)*csc(e + f*x)/(8*a**2*f*(a + b*sec(e + f*x)**2 - b)) - b*(9*a - 12*b)*sec(e + f*x)/(8*a**3*f*(a + b*sec(e + f*x)**2 - b)) + sqrt(b)*(-3*a + 6*b)*sqrt(a - b)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(2*a**4*f) - (3*a**2 - 24*a*b + 24*b**2)*atanh(cos(e + f*x))/(8*a**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_74():
    f = sin(e + f*x)**4/(a + b*tan(e + f*x)**2)**2
    F = -3*sqrt(a)*sqrt(b)*(a + b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(2*f*(a - b)**4) - 3*b*(3*a + b)*tan(e + f*x)/(8*f*(a - b)**3*(a + b*tan(e + f*x)**2)) + x*(3*a**2 + 18*a*b + 3*b**2)/(8*(a - b)**4) + sin(e + f*x)*cos(e + f*x)**3/(f*(a + b*tan(e + f*x)**2)*(4*a - 4*b)) - (5*a + b)*sin(e + f*x)*cos(e + f*x)/(8*f*(a - b)**2*(a + b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_75():
    f = sin(e + f*x)**2/(a + b*tan(e + f*x)**2)**2
    F = -b*tan(e + f*x)/(f*(a - b)**2*(a + b*tan(e + f*x)**2)) + x*(a + 3*b)/(2*(a - b)**3) - sin(e + f*x)*cos(e + f*x)/(f*(a + b*tan(e + f*x)**2)*(2*a - 2*b)) - sqrt(b)*(3*a + b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(2*sqrt(a)*f*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_76():
    f = (a + b*tan(e + f*x)**2)**(-2)
    F = x/(a - b)**2 - b*tan(e + f*x)/(2*a*f*(a - b)*(a + b*tan(e + f*x)**2)) - sqrt(b)*(3*a - b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_77():
    f = csc(e + f*x)**2/(a + b*tan(e + f*x)**2)**2
    F = cot(e + f*x)/(2*a*f*(a + b*tan(e + f*x)**2)) - 3*cot(e + f*x)/(2*a**2*f) - 3*sqrt(b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(2*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_78():
    f = csc(e + f*x)**4/(a + b*tan(e + f*x)**2)**2
    F = -cot(e + f*x)**3/(3*a**2*f) - b*(a - b)*tan(e + f*x)/(2*a**3*f*(a + b*tan(e + f*x)**2)) - (a - 2*b)*cot(e + f*x)/(a**3*f) - sqrt(b)*(3*a - 5*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(2*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_79():
    f = csc(e + f*x)**6/(a + b*tan(e + f*x)**2)**2
    F = -cot(e + f*x)**5/(5*a*f*(a + b*tan(e + f*x)**2)) - (10*a - 7*b)*cot(e + f*x)**3/(15*a**3*f) - b*(5*a**2 - 10*a*b + 7*b**2)*tan(e + f*x)/(10*a**4*f*(a + b*tan(e + f*x)**2)) - (5*a**2 - 20*a*b + 14*b**2)*cot(e + f*x)/(5*a**4*f) - sqrt(b)*(a - b)*(3*a - 7*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(2*a**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_80():
    f = sin(e + f*x)**5/(a + b*tan(e + f*x)**2)**3
    F = -sqrt(b)*(15*a**2 + 40*a*b + 8*b**2)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(8*f*(a - b)**(sympy.S(11)/2)) - b*(5*a**2 + 4*b**2)*sec(e + f*x)/(20*f*(a - b)**4*(a + b*sec(e + f*x)**2 - b)**2) - b*(35*a**2 + 40*a*b + 24*b**2)*sec(e + f*x)/(40*f*(a - b)**5*(a + b*sec(e + f*x)**2 - b)) - cos(e + f*x)**5/(f*(5*a - 5*b)*(a + b*sec(e + f*x)**2 - b)**2) + (10*a - b)*cos(e + f*x)**3/(15*f*(a - b)**4) - (5*a**2 + 20*a*b + 2*b**2)*cos(e + f*x)/(5*f*(a - b)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_81():
    f = sin(e + f*x)**3/(a + b*tan(e + f*x)**2)**3
    F = -a*b*sec(e + f*x)/(4*f*(a - b)**3*(a + b*sec(e + f*x)**2 - b)**2) - 5*sqrt(b)*(3*a + 4*b)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(8*f*(a - b)**(sympy.S(9)/2)) - b*(7*a + 4*b)*sec(e + f*x)/(8*f*(a - b)**4*(a + b*sec(e + f*x)**2 - b)) + cos(e + f*x)**3/(3*f*(a - b)**3) - (a + 2*b)*cos(e + f*x)/(f*(a - b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_82():
    f = sin(e + f*x)/(a + b*tan(e + f*x)**2)**3
    F = -15*sqrt(b)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(8*f*(a - b)**(sympy.S(7)/2)) + cos(e + f*x)/(f*(4*a - 4*b)*(a + b*sec(e + f*x)**2 - b)**2) + 5*cos(e + f*x)/(8*f*(a - b)**2*(a + b*sec(e + f*x)**2 - b)) - 15*cos(e + f*x)/(8*f*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_83():
    f = csc(e + f*x)/(a + b*tan(e + f*x)**2)**3
    F = -b*sec(e + f*x)/(4*a*f*(a - b)*(a + b*sec(e + f*x)**2 - b)**2) - b*(7*a - 4*b)*sec(e + f*x)/(8*a**2*f*(a - b)**2*(a + b*sec(e + f*x)**2 - b)) - sqrt(b)*(15*a**2 - 20*a*b + 8*b**2)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(8*a**3*f*(a - b)**(sympy.S(5)/2)) - atanh(cos(e + f*x))/(a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_84():
    f = csc(e + f*x)**3/(a + b*tan(e + f*x)**2)**3
    F = -cot(e + f*x)*csc(e + f*x)/(2*a*f*(a + b*sec(e + f*x)**2 - b)**2) - 3*b*sec(e + f*x)/(4*a**2*f*(a + b*sec(e + f*x)**2 - b)**2) - b*(11*a - 12*b)*sec(e + f*x)/(8*a**3*f*(a - b)*(a + b*sec(e + f*x)**2 - b)) - sqrt(b)*(15*a**2 - 40*a*b + 24*b**2)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(8*a**4*f*(a - b)**(sympy.S(3)/2)) - (a - 6*b)*atanh(cos(e + f*x))/(2*a**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_85():
    f = csc(e + f*x)**5/(a + b*tan(e + f*x)**2)**3
    F = -cot(e + f*x)**3*csc(e + f*x)/(4*a*f*(a + b*sec(e + f*x)**2 - b)**2) - (5*a - 8*b)*cot(e + f*x)*csc(e + f*x)/(8*a**2*f*(a + b*sec(e + f*x)**2 - b)**2) - b*(7*a - 12*b)*sec(e + f*x)/(8*a**3*f*(a + b*sec(e + f*x)**2 - b)**2) - b*(3*a - 6*b)*sec(e + f*x)/(2*a**4*f*(a + b*sec(e + f*x)**2 - b)) - 3*sqrt(b)*(5*a**2 - 20*a*b + 16*b**2)*atan(sqrt(b)*sec(e + f*x)/sqrt(a - b))/(8*a**5*f*sqrt(a - b)) - (3*a**2 - 36*a*b + 48*b**2)*atanh(cos(e + f*x))/(8*a**5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_86():
    f = sin(e + f*x)**4/(a + b*tan(e + f*x)**2)**3
    F = -b*(7*a + 5*b)*tan(e + f*x)/(8*f*(a - b)**3*(a + b*tan(e + f*x)**2)**2) - 3*b*(a + b)*tan(e + f*x)/(2*f*(a - b)**4*(a + b*tan(e + f*x)**2)) + x*(3*a**2 + 30*a*b + 15*b**2)/(8*(a - b)**5) + sin(e + f*x)*cos(e + f*x)**3/(f*(a + b*tan(e + f*x)**2)**2*(4*a - 4*b)) - (5*a + 3*b)*sin(e + f*x)*cos(e + f*x)/(8*f*(a - b)**2*(a + b*tan(e + f*x)**2)**2) - 3*sqrt(b)*(5*a**2 + 10*a*b + b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(8*sqrt(a)*f*(a - b)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_87():
    f = sin(e + f*x)**2/(a + b*tan(e + f*x)**2)**3
    F = -3*b*tan(e + f*x)/(4*f*(a - b)**2*(a + b*tan(e + f*x)**2)**2) + x*(a + 5*b)/(2*(a - b)**4) - sin(e + f*x)*cos(e + f*x)/(f*(a + b*tan(e + f*x)**2)**2*(2*a - 2*b)) - b*(11*a + b)*tan(e + f*x)/(8*a*f*(a - b)**3*(a + b*tan(e + f*x)**2)) - sqrt(b)*(15*a**2 + 10*a*b - b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(8*a**(sympy.S(3)/2)*f*(a - b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_88():
    f = (a + b*tan(e + f*x)**2)**(-3)
    F = x/(a - b)**3 - b*tan(e + f*x)/(4*a*f*(a - b)*(a + b*tan(e + f*x)**2)**2) - b*(7*a - 3*b)*tan(e + f*x)/(8*a**2*f*(a - b)**2*(a + b*tan(e + f*x)**2)) - sqrt(b)*(15*a**2 - 10*a*b + 3*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*f*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_89():
    f = csc(e + f*x)**2/(a + b*tan(e + f*x)**2)**3
    F = cot(e + f*x)/(4*a*f*(a + b*tan(e + f*x)**2)**2) + 5*cot(e + f*x)/(8*a**2*f*(a + b*tan(e + f*x)**2)) - 15*cot(e + f*x)/(8*a**3*f) - 15*sqrt(b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(8*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_90():
    f = csc(e + f*x)**4/(a + b*tan(e + f*x)**2)**3
    F = -b*(a - b)*tan(e + f*x)/(4*a**3*f*(a + b*tan(e + f*x)**2)**2) - cot(e + f*x)**3/(3*a**3*f) - b*(7*a - 11*b)*tan(e + f*x)/(8*a**4*f*(a + b*tan(e + f*x)**2)) - (a - 3*b)*cot(e + f*x)/(a**4*f) + sqrt(b)*(-15*a + 35*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(8*a**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_91():
    f = csc(e + f*x)**6/(a + b*tan(e + f*x)**2)**3
    F = -cot(e + f*x)**5/(5*a*f*(a + b*tan(e + f*x)**2)**2) - b*(5*a**2 - 10*a*b + 9*b**2)*tan(e + f*x)/(20*a**4*f*(a + b*tan(e + f*x)**2)**2) - (10*a - 9*b)*cot(e + f*x)**3/(15*a**4*f) - b*(35*a**2 - 110*a*b + 99*b**2)*tan(e + f*x)/(40*a**5*f*(a + b*tan(e + f*x)**2)) - (5*a**2 - 30*a*b + 27*b**2)*cot(e + f*x)/(5*a**5*f) - sqrt(b)*(15*a**2 - 70*a*b + 63*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(8*a**(sympy.S(11)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_92():
    f = sqrt(a + b*tan(e + f*x)**2)*sin(e + f*x)**5
    F = sqrt(b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/f - sqrt(a + b*sec(e + f*x)**2 - b)*cos(e + f*x)/f - (a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)*cos(e + f*x)**5/(f*(5*a - 5*b)) + (10*a - 8*b)*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)*cos(e + f*x)**3/(15*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_93():
    f = sqrt(a + b*tan(e + f*x)**2)*sin(e + f*x)**3
    F = sqrt(b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/f - sqrt(a + b*sec(e + f*x)**2 - b)*cos(e + f*x)/f + (a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)*cos(e + f*x)**3/(f*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_94():
    f = sqrt(a + b*tan(e + f*x)**2)*sin(e + f*x)
    F = sqrt(b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/f - sqrt(a + b*sec(e + f*x)**2 - b)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_95():
    f = sqrt(a + b*tan(e + f*x)**2)*csc(e + f*x)
    F = -sqrt(a)*atanh(sqrt(a)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/f + sqrt(b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_96():
    f = sqrt(a + b*tan(e + f*x)**2)*csc(e + f*x)**3
    F = sqrt(b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/f - sqrt(a + b*sec(e + f*x)**2 - b)*cot(e + f*x)*csc(e + f*x)/(2*f) - (a + b)*atanh(sqrt(a)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(2*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_97():
    f = sqrt(a + b*tan(e + f*x)**2)*csc(e + f*x)**5
    F = sqrt(b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/f - sqrt(a + b*sec(e + f*x)**2 - b)*cot(e + f*x)*csc(e + f*x)**3/(4*f) - (3*a + b)*sqrt(a + b*sec(e + f*x)**2 - b)*cot(e + f*x)*csc(e + f*x)/(8*a*f) - (3*a**2 + 6*a*b - b**2)*atanh(sqrt(a)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(8*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_98():
    f = sqrt(a + b*tan(e + f*x)**2)*sin(e + f*x)**4
    F = sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f - sqrt(a + b*tan(e + f*x)**2)*(3*a - 4*b)*sin(e + f*x)*cos(e + f*x)/(f*(8*a - 8*b)) - sqrt(a + b*tan(e + f*x)**2)*sin(e + f*x)**3*cos(e + f*x)/(4*f) + (3*a**2 - 12*a*b + 8*b**2)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(8*f*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_99():
    f = sqrt(a + b*tan(e + f*x)**2)*sin(e + f*x)**2
    F = sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f + (a - 2*b)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(2*f*sqrt(a - b)) - sqrt(a + b*tan(e + f*x)**2)*sin(e + f*x)*cos(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_100():
    f = sqrt(a + b*tan(e + f*x)**2)
    F = sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f + sqrt(a - b)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_101():
    f = sqrt(a + b*tan(e + f*x)**2)*csc(e + f*x)**2
    F = sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f - sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_102():
    f = sqrt(a + b*tan(e + f*x)**2)*csc(e + f*x)**4
    F = sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f - sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)/f - (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**3/(3*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_103():
    f = sqrt(a + b*tan(e + f*x)**2)*csc(e + f*x)**6
    F = sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f - sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)/f - (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**5/(5*a*f) - (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*(10*a - 2*b)*cot(e + f*x)**3/(15*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_104():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)**5
    F = sqrt(b)*(3*a - 7*b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(2*f) + b*(3*a - 7*b)*sqrt(a + b*sec(e + f*x)**2 - b)*sec(e + f*x)/(f*(2*a - 2*b)) - (3*a - 7*b)*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)*cos(e + f*x)/(f*(3*a - 3*b)) - (a + b*sec(e + f*x)**2 - b)**(sympy.S(5)/2)*cos(e + f*x)**5/(f*(5*a - 5*b)) + 2*(a + b*sec(e + f*x)**2 - b)**(sympy.S(5)/2)*cos(e + f*x)**3/(f*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_105():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)**3
    F = sqrt(b)*(3*a - 5*b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(2*f) + b*(3*a - 5*b)*sqrt(a + b*sec(e + f*x)**2 - b)*sec(e + f*x)/(f*(2*a - 2*b)) - (3*a - 5*b)*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)*cos(e + f*x)/(f*(3*a - 3*b)) + (a + b*sec(e + f*x)**2 - b)**(sympy.S(5)/2)*cos(e + f*x)**3/(f*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_106():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)
    F = sqrt(b)*(3*a - 3*b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(2*f) + 3*b*sqrt(a + b*sec(e + f*x)**2 - b)*sec(e + f*x)/(2*f) - (a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_107():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/f + sqrt(b)*(3*a - b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(2*f) + b*sqrt(a + b*sec(e + f*x)**2 - b)*sec(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_108():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**3
    F = -sqrt(a)*(a + 3*b)*atanh(sqrt(a)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(2*f) + sqrt(b)*(3*a + b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(2*f) + b*sqrt(a + b*sec(e + f*x)**2 - b)*sec(e + f*x)/f - (a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)*cot(e + f*x)*csc(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_109():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**5
    F = 3*sqrt(b)*(a + b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(2*f) - (3*a + 3*b)*sqrt(a + b*sec(e + f*x)**2 - b)*csc(e + f*x)**2*sec(e + f*x)/(8*f) + (3*a + 9*b)*sqrt(a + b*sec(e + f*x)**2 - b)*sec(e + f*x)/(8*f) - (a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)*cot(e + f*x)*csc(e + f*x)**3/(4*f) + (-3*a**2 - 18*a*b - 3*b**2)*atanh(sqrt(a)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(8*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_110():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)**4
    F = sqrt(b)*(3*a - 6*b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(2*f) - (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)**3*cos(e + f*x)/(4*f) - sqrt(a + b*tan(e + f*x)**2)*(3*a - 12*b)*tan(e + f*x)/(8*f) + sqrt(a + b*tan(e + f*x)**2)*(3*a - 6*b)*sin(e + f*x)**2*tan(e + f*x)/(8*f) + (3*a**2 - 24*a*b + 24*b**2)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(8*f*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_111():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)**2
    F = sqrt(b)*(3*a - 4*b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(2*f) + b*sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)/f + (a - 4*b)*sqrt(a - b)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(2*f) - (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)*cos(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_112():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = sqrt(b)*(3*a - 2*b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(2*f) + b*sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)/(2*f) + (a - b)**(sympy.S(3)/2)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_113():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**2
    F = 3*a*sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(2*f) + 3*b*sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)/(2*f) - (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_114():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**4
    F = sqrt(b)*(3*a + 2*b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(2*f) + b*sqrt(a + b*tan(e + f*x)**2)*(3*a + 2*b)*tan(e + f*x)/(2*a*f) - (a + b*tan(e + f*x)**2)**(sympy.S(5)/2)*cot(e + f*x)**3/(3*a*f) - (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*(3*a + 2*b)*cot(e + f*x)/(3*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_115():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**6
    F = sqrt(b)*(3*a + 4*b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(2*f) + b*sqrt(a + b*tan(e + f*x)**2)*(3*a + 4*b)*tan(e + f*x)/(2*a*f) - (a + b*tan(e + f*x)**2)**(sympy.S(5)/2)*cot(e + f*x)**5/(5*a*f) - 2*(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)*cot(e + f*x)**3/(3*a*f) - (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*(3*a + 4*b)*cot(e + f*x)/(3*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_116():
    f = sin(e + f*x)**5/sqrt(a + b*tan(e + f*x)**2)
    F = -sqrt(a + b*sec(e + f*x)**2 - b)*cos(e + f*x)**5/(f*(5*a - 5*b)) + (10*a - 6*b)*sqrt(a + b*sec(e + f*x)**2 - b)*cos(e + f*x)**3/(15*f*(a - b)**2) - sqrt(a + b*sec(e + f*x)**2 - b)*(15*a**2 - 10*a*b + 3*b**2)*cos(e + f*x)/(15*f*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_117():
    f = sin(e + f*x)**3/sqrt(a + b*tan(e + f*x)**2)
    F = sqrt(a + b*sec(e + f*x)**2 - b)*cos(e + f*x)**3/(f*(3*a - 3*b)) - (3*a - b)*sqrt(a + b*sec(e + f*x)**2 - b)*cos(e + f*x)/(3*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_118():
    f = sin(e + f*x)/sqrt(a + b*tan(e + f*x)**2)
    F = -sqrt(a + b*sec(e + f*x)**2 - b)*cos(e + f*x)/(f*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_119():
    f = csc(e + f*x)/sqrt(a + b*tan(e + f*x)**2)
    F = -atanh(sqrt(a)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_120():
    f = csc(e + f*x)**3/sqrt(a + b*tan(e + f*x)**2)
    F = -sqrt(a + b*sec(e + f*x)**2 - b)*cot(e + f*x)*csc(e + f*x)/(2*a*f) - (a - b)*atanh(sqrt(a)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(2*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_121():
    f = csc(e + f*x)**5/sqrt(a + b*tan(e + f*x)**2)
    F = -sqrt(a + b*sec(e + f*x)**2 - b)*cot(e + f*x)**3*csc(e + f*x)/(4*a*f) - (5*a - 3*b)*sqrt(a + b*sec(e + f*x)**2 - b)*cot(e + f*x)*csc(e + f*x)/(8*a**2*f) - 3*(a - b)**2*atanh(sqrt(a)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(8*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_122():
    f = sin(e + f*x)**4/sqrt(a + b*tan(e + f*x)**2)
    F = 3*a**2*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(8*f*(a - b)**(sympy.S(5)/2)) + sqrt(a + b*tan(e + f*x)**2)*sin(e + f*x)*cos(e + f*x)**3/(f*(4*a - 4*b)) - sqrt(a + b*tan(e + f*x)**2)*(5*a - 2*b)*sin(e + f*x)*cos(e + f*x)/(8*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_123():
    f = sin(e + f*x)**2/sqrt(a + b*tan(e + f*x)**2)
    F = a*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(2*f*(a - b)**(sympy.S(3)/2)) - sqrt(a + b*tan(e + f*x)**2)*sin(e + f*x)*cos(e + f*x)/(f*(2*a - 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_124():
    f = 1/sqrt(a + b*tan(e + f*x)**2)
    F = atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_125():
    f = csc(e + f*x)**2/sqrt(a + b*tan(e + f*x)**2)
    F = -sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_126():
    f = csc(e + f*x)**4/sqrt(a + b*tan(e + f*x)**2)
    F = -sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**3/(3*a*f) - sqrt(a + b*tan(e + f*x)**2)*(3*a - 2*b)*cot(e + f*x)/(3*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_127():
    f = csc(e + f*x)**6/sqrt(a + b*tan(e + f*x)**2)
    F = -sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**5/(5*a*f) - sqrt(a + b*tan(e + f*x)**2)*(10*a - 4*b)*cot(e + f*x)**3/(15*a**2*f) - sqrt(a + b*tan(e + f*x)**2)*(15*a**2 - 20*a*b + 8*b**2)*cot(e + f*x)/(15*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_128():
    f = sin(e + f*x)**5/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -2*b*(15*a**2 + 10*a*b - b**2)*sec(e + f*x)/(15*f*(a - b)**4*sqrt(a + b*sec(e + f*x)**2 - b)) - cos(e + f*x)**5/(f*(5*a - 5*b)*sqrt(a + b*sec(e + f*x)**2 - b)) + (10*a - 4*b)*cos(e + f*x)**3/(15*f*(a - b)**2*sqrt(a + b*sec(e + f*x)**2 - b)) - (15*a**2 + 10*a*b - b**2)*cos(e + f*x)/(15*f*(a - b)**3*sqrt(a + b*sec(e + f*x)**2 - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_129():
    f = sin(e + f*x)**3/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -2*b*(3*a + b)*sec(e + f*x)/(3*f*(a - b)**3*sqrt(a + b*sec(e + f*x)**2 - b)) + cos(e + f*x)**3/(f*(3*a - 3*b)*sqrt(a + b*sec(e + f*x)**2 - b)) - (3*a + b)*cos(e + f*x)/(3*f*(a - b)**2*sqrt(a + b*sec(e + f*x)**2 - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_130():
    f = sin(e + f*x)/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -2*b*sec(e + f*x)/(f*(a - b)**2*sqrt(a + b*sec(e + f*x)**2 - b)) - cos(e + f*x)/(f*(a - b)*sqrt(a + b*sec(e + f*x)**2 - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_131():
    f = csc(e + f*x)/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -b*sec(e + f*x)/(a*f*(a - b)*sqrt(a + b*sec(e + f*x)**2 - b)) - atanh(sqrt(a)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_132():
    f = csc(e + f*x)**3/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -cot(e + f*x)*csc(e + f*x)/(2*a*f*sqrt(a + b*sec(e + f*x)**2 - b)) - 3*b*sec(e + f*x)/(2*a**2*f*sqrt(a + b*sec(e + f*x)**2 - b)) - (a - 3*b)*atanh(sqrt(a)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(2*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_133():
    f = csc(e + f*x)**5/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -cot(e + f*x)**3*csc(e + f*x)/(4*a*f*sqrt(a + b*sec(e + f*x)**2 - b)) - (5*a - 5*b)*cot(e + f*x)*csc(e + f*x)/(8*a**2*f*sqrt(a + b*sec(e + f*x)**2 - b)) - b*(13*a - 15*b)*sec(e + f*x)/(8*a**3*f*sqrt(a + b*sec(e + f*x)**2 - b)) + (-3*a + 15*b)*(a - b)*atanh(sqrt(a)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(8*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_134():
    f = sin(e + f*x)**4/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -5*a*sin(e + f*x)*cos(e + f*x)/(8*f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2)) + 3*a*(a + 4*b)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(8*f*(a - b)**(sympy.S(7)/2)) - b*(13*a + 2*b)*tan(e + f*x)/(8*f*(a - b)**3*sqrt(a + b*tan(e + f*x)**2)) + sin(e + f*x)*cos(e + f*x)**3/(f*sqrt(a + b*tan(e + f*x)**2)*(4*a - 4*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_135():
    f = sin(e + f*x)**2/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -3*b*tan(e + f*x)/(2*f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2)) - sin(e + f*x)*cos(e + f*x)/(f*sqrt(a + b*tan(e + f*x)**2)*(2*a - 2*b)) + (a + 2*b)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(2*f*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_136():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(-3)/2)
    F = atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(3)/2)) - b*tan(e + f*x)/(a*f*(a - b)*sqrt(a + b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_137():
    f = csc(e + f*x)**2/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -cot(e + f*x)/(a*f*sqrt(a + b*tan(e + f*x)**2)) - 2*b*tan(e + f*x)/(a**2*f*sqrt(a + b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_138():
    f = csc(e + f*x)**4/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -cot(e + f*x)**3/(3*a*f*sqrt(a + b*tan(e + f*x)**2)) - (3*a - 4*b)*cot(e + f*x)/(3*a**2*f*sqrt(a + b*tan(e + f*x)**2)) - b*(6*a - 8*b)*tan(e + f*x)/(3*a**3*f*sqrt(a + b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_139():
    f = csc(e + f*x)**6/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -cot(e + f*x)**5/(5*a*f*sqrt(a + b*tan(e + f*x)**2)) - (10*a - 6*b)*cot(e + f*x)**3/(15*a**2*f*sqrt(a + b*tan(e + f*x)**2)) - (15*a**2 - 40*a*b + 24*b**2)*cot(e + f*x)/(15*a**3*f*sqrt(a + b*tan(e + f*x)**2)) - 2*b*(15*a**2 - 40*a*b + 24*b**2)*tan(e + f*x)/(15*a**4*f*sqrt(a + b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_140():
    f = sin(e + f*x)**5/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -4*b*(5*a**2 + 10*a*b + b**2)*sec(e + f*x)/(15*f*(a - b)**4*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)) - 8*b*(5*a**2 + 10*a*b + b**2)*sec(e + f*x)/(15*f*(a - b)**5*sqrt(a + b*sec(e + f*x)**2 - b)) - cos(e + f*x)**5/(f*(5*a - 5*b)*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)) + (10*a - 2*b)*cos(e + f*x)**3/(15*f*(a - b)**2*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)) - (5*a**2 + 10*a*b + b**2)*cos(e + f*x)/(5*f*(a - b)**3*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_141():
    f = sin(e + f*x)**3/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -4*b*(a + b)*sec(e + f*x)/(3*f*(a - b)**3*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)) - 8*b*(a + b)*sec(e + f*x)/(3*f*(a - b)**4*sqrt(a + b*sec(e + f*x)**2 - b)) + cos(e + f*x)**3/(f*(3*a - 3*b)*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)) - (a + b)*cos(e + f*x)/(f*(a - b)**2*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_142():
    f = sin(e + f*x)/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -4*b*sec(e + f*x)/(3*f*(a - b)**2*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)) - 8*b*sec(e + f*x)/(3*f*(a - b)**3*sqrt(a + b*sec(e + f*x)**2 - b)) - cos(e + f*x)/(f*(a - b)*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_143():
    f = csc(e + f*x)/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -b*sec(e + f*x)/(3*a*f*(a - b)*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)) - b*(5*a - 3*b)*sec(e + f*x)/(3*a**2*f*(a - b)**2*sqrt(a + b*sec(e + f*x)**2 - b)) - atanh(sqrt(a)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_144():
    f = csc(e + f*x)**3/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -cot(e + f*x)*csc(e + f*x)/(2*a*f*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)) - 5*b*sec(e + f*x)/(6*a**2*f*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)) - b*(13*a - 15*b)*sec(e + f*x)/(6*a**3*f*(a - b)*sqrt(a + b*sec(e + f*x)**2 - b)) - (a - 5*b)*atanh(sqrt(a)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(2*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_145():
    f = csc(e + f*x)**5/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -cot(e + f*x)**3*csc(e + f*x)/(4*a*f*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)) - (5*a - 7*b)*cot(e + f*x)*csc(e + f*x)/(8*a**2*f*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)) - b*(23*a - 35*b)*sec(e + f*x)/(24*a**3*f*(a + b*sec(e + f*x)**2 - b)**(sympy.S(3)/2)) - b*(55*a - 105*b)*sec(e + f*x)/(24*a**4*f*sqrt(a + b*sec(e + f*x)**2 - b)) - (3*a**2 - 30*a*b + 35*b**2)*atanh(sqrt(a)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2 - b))/(8*a**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_146():
    f = sin(e + f*x)**4/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -b*(23*a + 12*b)*tan(e + f*x)/(24*f*(a - b)**3*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - 5*b*(11*a + 10*b)*tan(e + f*x)/(24*f*(a - b)**4*sqrt(a + b*tan(e + f*x)**2)) + sin(e + f*x)*cos(e + f*x)**3/(f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*(4*a - 4*b)) - (5*a + 2*b)*sin(e + f*x)*cos(e + f*x)/(8*f*(a - b)**2*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) + (3*a**2 + 24*a*b + 8*b**2)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(8*f*(a - b)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_147():
    f = sin(e + f*x)**2/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -5*b*tan(e + f*x)/(6*f*(a - b)**2*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - sin(e + f*x)*cos(e + f*x)/(f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*(2*a - 2*b)) + (a + 4*b)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(2*f*(a - b)**(sympy.S(7)/2)) - b*(13*a + 2*b)*tan(e + f*x)/(6*a*f*(a - b)**3*sqrt(a + b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_148():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(-5)/2)
    F = atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(5)/2)) - b*tan(e + f*x)/(3*a*f*(a - b)*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - b*(5*a - 2*b)*tan(e + f*x)/(3*a**2*f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_149():
    f = csc(e + f*x)**2/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -cot(e + f*x)/(a*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - 4*b*tan(e + f*x)/(3*a**2*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - 8*b*tan(e + f*x)/(3*a**3*f*sqrt(a + b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_150():
    f = csc(e + f*x)**4/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -cot(e + f*x)**3/(3*a*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - (a - 2*b)*cot(e + f*x)/(a**2*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - b*(4*a - 8*b)*tan(e + f*x)/(3*a**3*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - b*(8*a - 16*b)*tan(e + f*x)/(3*a**4*f*sqrt(a + b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_151():
    f = csc(e + f*x)**6/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -cot(e + f*x)**5/(5*a*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - (10*a - 8*b)*cot(e + f*x)**3/(15*a**2*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - (5*a**2 - 20*a*b + 16*b**2)*cot(e + f*x)/(5*a**3*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - 4*b*(5*a**2 - 20*a*b + 16*b**2)*tan(e + f*x)/(15*a**4*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - 8*b*(5*a**2 - 20*a*b + 16*b**2)*tan(e + f*x)/(15*a**5*f*sqrt(a + b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_152():
    f = (b*tan(e + f*x)**2)**p*(d*sin(e + f*x))**m
    F = (b*tan(e + f*x)**2)**p*(d*sin(e + f*x))**m*(cos(e + f*x)**2)**(p + sympy.S.Half)*tan(e + f*x)*hyper((p + sympy.S.Half, m/2 + p + sympy.S.Half), (m/2 + p + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(m + 2*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_153():
    f = (d*sin(e + f*x))**m*(a + b*tan(e + f*x)**2)**p
    F = (d*sin(e + f*x))**m*(a + b*tan(e + f*x)**2)**p*(sec(e + f*x)**2)**(m/2)*tan(e + f*x)*appellf1(m/2 + sympy.S.Half, -p, m/2 + 1, m/2 + sympy.S(3)/2, -b*tan(e + f*x)**2/a, -tan(e + f*x)**2)/(f*(1 + b*tan(e + f*x)**2/a)**p*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_154():
    f = (a + b*tan(e + f*x)**2)**p*sin(e + f*x)**5
    F = -(a + b*sec(e + f*x)**2 - b)**(p + 1)*cos(e + f*x)**5/(f*(5*a - 5*b)) + (a + b*sec(e + f*x)**2 - b)**(p + 1)*(10*a - 2*b*p - 7*b)*cos(e + f*x)**3/(15*f*(a - b)**2) - (a + b*sec(e + f*x)**2 - b)**p*(15*a**2 - 20*a*b*(p + 1) + 4*b**2*(p**2 + 3*p + 2))*cos(e + f*x)*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*sec(e + f*x)**2/(a - b))/(15*f*(a - b)**2*(b*sec(e + f*x)**2/(a - b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_155():
    f = (a + b*tan(e + f*x)**2)**p*sin(e + f*x)**3
    F = -(3*a - 2*b*(p + 1))*(a + b*sec(e + f*x)**2 - b)**p*cos(e + f*x)*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*sec(e + f*x)**2/(a - b))/(f*(3*a - 3*b)*(b*sec(e + f*x)**2/(a - b) + 1)**p) + (a + b*sec(e + f*x)**2 - b)**(p + 1)*cos(e + f*x)**3/(f*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_156():
    f = (a + b*tan(e + f*x)**2)**p*sin(e + f*x)
    F = -(a + b*sec(e + f*x)**2 - b)**p*cos(e + f*x)*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*sec(e + f*x)**2/(a - b))/(f*(b*sec(e + f*x)**2/(a - b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_157():
    f = (a + b*tan(e + f*x)**2)**p*csc(e + f*x)
    F = -(a + b*sec(e + f*x)**2 - b)**p*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, sec(e + f*x)**2, -b*sec(e + f*x)**2/(a - b))*sec(e + f*x)/(f*(b*sec(e + f*x)**2/(a - b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_158():
    f = (a + b*tan(e + f*x)**2)**p*csc(e + f*x)**3
    F = (a + b*sec(e + f*x)**2 - b)**p*appellf1(sympy.S(3)/2, 2, -p, sympy.S(5)/2, sec(e + f*x)**2, -b*sec(e + f*x)**2/(a - b))*sec(e + f*x)**3/(3*f*(b*sec(e + f*x)**2/(a - b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_159():
    f = (a + b*tan(e + f*x)**2)**p*sin(e + f*x)**2
    F = (a + b*tan(e + f*x)**2)**p*tan(e + f*x)**3*appellf1(sympy.S(3)/2, 2, -p, sympy.S(5)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/a)/(3*f*(1 + b*tan(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_160():
    f = (a + b*tan(e + f*x)**2)**p
    F = (a + b*tan(e + f*x)**2)**p*tan(e + f*x)*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/a)/(f*(1 + b*tan(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_161():
    f = (a + b*tan(e + f*x)**2)**p*csc(e + f*x)**2
    F = -(a + b*tan(e + f*x)**2)**p*cot(e + f*x)*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*tan(e + f*x)**2/a)/(f*(1 + b*tan(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_162():
    f = (a + b*tan(e + f*x)**2)**p*csc(e + f*x)**4
    F = -(a + b*tan(e + f*x)**2)**(p + 1)*cot(e + f*x)**3/(3*a*f) - (a + b*tan(e + f*x)**2)**p*(3*a - b*(1 - 2*p))*cot(e + f*x)*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*tan(e + f*x)**2/a)/(3*a*f*(1 + b*tan(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_163():
    f = (a + b*tan(e + f*x)**2)**p*csc(e + f*x)**6
    F = -(a + b*tan(e + f*x)**2)**(p + 1)*cot(e + f*x)**5/(5*a*f) - (a + b*tan(e + f*x)**2)**(p + 1)*(10*a - b*(3 - 2*p))*cot(e + f*x)**3/(15*a**2*f) - (a + b*tan(e + f*x)**2)**p*(15*a**2 - b*(1 - 2*p)*(10*a - b*(3 - 2*p)))*cot(e + f*x)*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*tan(e + f*x)**2/a)/(15*a**2*f*(1 + b*tan(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_164():
    f = (b*(c*tan(e + f*x))**n)**p*(d*sin(e + f*x))**m
    F = (b*(c*tan(e + f*x))**n)**p*(d*sin(e + f*x))**m*(cos(e + f*x)**2)**(n*p/2 + sympy.S.Half)*tan(e + f*x)*hyper((n*p/2 + sympy.S.Half, m/2 + n*p/2 + sympy.S.Half), (m/2 + n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(m + n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_165():
    f = (b*(c*tan(e + f*x))**n)**p*sin(e + f*x)**2
    F = (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)**3*hyper((2, n*p/2 + sympy.S(3)/2), (n*p/2 + sympy.S(5)/2,), -tan(e + f*x)**2)/(f*(n*p + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_166():
    f = (b*(c*tan(e + f*x))**n)**p
    F = (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)*hyper((1, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_167():
    f = (b*(c*tan(e + f*x))**n)**p*csc(e + f*x)**2
    F = -(b*(c*tan(e + f*x))**n)**p*cot(e + f*x)/(f*(-n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_168():
    f = (b*(c*tan(e + f*x))**n)**p*csc(e + f*x)**4
    F = -(b*(c*tan(e + f*x))**n)**p*cot(e + f*x)**3/(f*(-n*p + 3)) - (b*(c*tan(e + f*x))**n)**p*cot(e + f*x)/(f*(-n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_169():
    f = (b*(c*tan(e + f*x))**n)**p*csc(e + f*x)**6
    F = -(b*(c*tan(e + f*x))**n)**p*cot(e + f*x)**5/(f*(-n*p + 5)) - 2*(b*(c*tan(e + f*x))**n)**p*cot(e + f*x)**3/(f*(-n*p + 3)) - (b*(c*tan(e + f*x))**n)**p*cot(e + f*x)/(f*(-n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_170():
    f = (b*(c*tan(e + f*x))**n)**p*sin(e + f*x)**3
    F = (b*(c*tan(e + f*x))**n)**p*(cos(e + f*x)**2)**(n*p/2 + sympy.S.Half)*sin(e + f*x)**3*tan(e + f*x)*hyper((n*p/2 + sympy.S.Half, n*p/2 + 2), (n*p/2 + 3,), sin(e + f*x)**2)/(f*(n*p + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_171():
    f = (b*(c*tan(e + f*x))**n)**p*sin(e + f*x)
    F = (b*(c*tan(e + f*x))**n)**p*(cos(e + f*x)**2)**(n*p/2 + sympy.S.Half)*sin(e + f*x)*tan(e + f*x)*hyper((n*p/2 + sympy.S.Half, n*p/2 + 1), (n*p/2 + 2,), sin(e + f*x)**2)/(f*(n*p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_172():
    f = (b*(c*tan(e + f*x))**n)**p*csc(e + f*x)
    F = (b*(c*tan(e + f*x))**n)**p*(cos(e + f*x)**2)**(n*p/2 + sympy.S.Half)*hyper((n*p/2, n*p/2 + sympy.S.Half), (n*p/2 + 1,), sin(e + f*x)**2)*sec(e + f*x)/(f*n*p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_173():
    f = (b*(c*tan(e + f*x))**n)**p*csc(e + f*x)**3
    F = -(b*(c*tan(e + f*x))**n)**p*(cos(e + f*x)**2)**(n*p/2 + sympy.S.Half)*csc(e + f*x)**2*hyper((n*p/2 - 1, n*p/2 + sympy.S.Half), (n*p/2,), sin(e + f*x)**2)*sec(e + f*x)/(f*(-n*p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_174():
    f = (d*sin(e + f*x))**m*(a + b*tan(e + f*x)**n)**p
    F = sympy.Function('Unintegrable')((((Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_175():
    f = (b*tan(e + f*x)**2)**p*(d*cos(e + f*x))**m
    F = (b*tan(e + f*x)**2)**p*(d*cos(e + f*x))**m*(cos(e + f*x)**2)**(-m/2 + p + sympy.S.Half)*tan(e + f*x)*hyper((p + sympy.S.Half, -m/2 + p + sympy.S.Half), (p + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_176():
    f = (d*cos(e + f*x))**m*(a + b*tan(e + f*x)**2)**p
    F = (d*cos(e + f*x))**m*(a + b*tan(e + f*x)**2)**p*(sec(e + f*x)**2)**(m/2)*tan(e + f*x)*appellf1(sympy.S.Half, -p, m/2 + 1, sympy.S(3)/2, -b*tan(e + f*x)**2/a, -tan(e + f*x)**2)/(f*(1 + b*tan(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_177():
    f = (b*(c*tan(e + f*x))**n)**p*(d*cos(e + f*x))**m
    F = (b*(c*tan(e + f*x))**n)**p*(d*cos(e + f*x))**m*(cos(e + f*x)**2)**(-m/2 + n*p/2 + sympy.S.Half)*tan(e + f*x)*hyper((n*p/2 + sympy.S.Half, -m/2 + n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_178():
    f = (d*cos(e + f*x))**m*(a + b*(c*tan(e + f*x))**n)**p
    F = ((Symbol('d') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')) * sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p')) * (((sympy.sec((Symbol('e') + (Symbol('f') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_179():
    f = (a*tan(c + d*x)**2 + a)**4
    F = a**4*tan(c + d*x)**7/(7*d) + 3*a**4*tan(c + d*x)**5/(5*d) + a**4*tan(c + d*x)**3/d + a**4*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_180():
    f = (a*tan(c + d*x)**2 + a)**3
    F = a**3*tan(c + d*x)**5/(5*d) + 2*a**3*tan(c + d*x)**3/(3*d) + a**3*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_181():
    f = (a*tan(c + d*x)**2 + a)**2
    F = a**2*tan(c + d*x)**3/(3*d) + a**2*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_182():
    f = 1/(a*tan(c + d*x)**2 + a)
    F = x/(2*a) + sin(c + d*x)*cos(c + d*x)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_183():
    f = (a*tan(c + d*x)**2 + a)**(-2)
    F = 3*x/(8*a**2) + sin(c + d*x)*cos(c + d*x)**3/(4*a**2*d) + 3*sin(c + d*x)*cos(c + d*x)/(8*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_184():
    f = (a*tan(c + d*x)**2 + a)**(-3)
    F = 5*x/(16*a**3) + sin(c + d*x)*cos(c + d*x)**5/(6*a**3*d) + 5*sin(c + d*x)*cos(c + d*x)**3/(24*a**3*d) + 5*sin(c + d*x)*cos(c + d*x)/(16*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_185():
    f = (a + b*tan(e + f*x)**2)*tan(e + f*x)**5
    F = b*tan(e + f*x)**6/(6*f) - (a - b)*log(cos(e + f*x))/f + (a - b)*tan(e + f*x)**4/(4*f) - (a - b)*tan(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_186():
    f = (a + b*tan(e + f*x)**2)*tan(e + f*x)**3
    F = b*tan(e + f*x)**4/(4*f) + (a - b)*log(cos(e + f*x))/f + (a - b)*tan(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_187():
    f = (a + b*tan(e + f*x)**2)*tan(e + f*x)
    F = b*tan(e + f*x)**2/(2*f) - (a - b)*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_188():
    f = (a + b*tan(e + f*x)**2)*cot(e + f*x)
    F = a*log(sin(e + f*x))/f - b*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_189():
    f = (a + b*tan(e + f*x)**2)*cot(e + f*x)**3
    F = -a*cot(e + f*x)**2/(2*f) - (a - b)*log(sin(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_190():
    f = (a + b*tan(e + f*x)**2)*cot(e + f*x)**5
    F = -a*cot(e + f*x)**4/(4*f) + (a - b)*log(sin(e + f*x))/f + (a - b)*cot(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_191():
    f = (a + b*tan(e + f*x)**2)*tan(e + f*x)**6
    F = b*tan(e + f*x)**7/(7*f) - x*(a - b) + (a - b)*tan(e + f*x)**5/(5*f) - (a - b)*tan(e + f*x)**3/(3*f) + (a - b)*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_192():
    f = (a + b*tan(e + f*x)**2)*tan(e + f*x)**4
    F = b*tan(e + f*x)**5/(5*f) + x*(a - b) + (a - b)*tan(e + f*x)**3/(3*f) - (a - b)*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_193():
    f = (a + b*tan(e + f*x)**2)*tan(e + f*x)**2
    F = b*tan(e + f*x)**3/(3*f) - x*(a - b) + (a - b)*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_194():
    f = a + b*tan(e + f*x)**2
    F = a*x - b*x + b*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_195():
    f = (a + b*tan(e + f*x)**2)*cot(e + f*x)**2
    F = -a*cot(e + f*x)/f - x*(a - b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_196():
    f = (a + b*tan(e + f*x)**2)*cot(e + f*x)**4
    F = -a*cot(e + f*x)**3/(3*f) + x*(a - b) + (a - b)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_197():
    f = (a + b*tan(e + f*x)**2)*cot(e + f*x)**6
    F = -a*cot(e + f*x)**5/(5*f) - x*(a - b) + (a - b)*cot(e + f*x)**3/(3*f) - (a - b)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_198():
    f = (a + b*tan(e + f*x)**2)**2*tan(e + f*x)**5
    F = b**2*tan(e + f*x)**8/(8*f) + b*(2*a - b)*tan(e + f*x)**6/(6*f) - (a - b)**2*log(cos(e + f*x))/f + (a - b)**2*tan(e + f*x)**4/(4*f) - (a - b)**2*tan(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_199():
    f = (a + b*tan(e + f*x)**2)**2*tan(e + f*x)**3
    F = b**2*tan(e + f*x)**6/(6*f) + b*(2*a - b)*tan(e + f*x)**4/(4*f) + (a - b)**2*log(cos(e + f*x))/f + (a - b)**2*tan(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_200():
    f = (a + b*tan(e + f*x)**2)**2*tan(e + f*x)
    F = b*(a - b)*tan(e + f*x)**2/(2*f) - (a - b)**2*log(cos(e + f*x))/f + (a + b*tan(e + f*x)**2)**2/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_201():
    f = (a + b*tan(e + f*x)**2)**2*cot(e + f*x)
    F = a**2*log(tan(e + f*x))/f + b**2*tan(e + f*x)**2/(2*f) + (a - b)**2*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_202():
    f = (a + b*tan(e + f*x)**2)**2*cot(e + f*x)**3
    F = -a**2*cot(e + f*x)**2/(2*f) - a*(a - 2*b)*log(tan(e + f*x))/f - (a - b)**2*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_203():
    f = (a + b*tan(e + f*x)**2)**2*cot(e + f*x)**5
    F = -a**2*cot(e + f*x)**4/(4*f) + a*(a - 2*b)*cot(e + f*x)**2/(2*f) + (a - b)**2*log(cos(e + f*x))/f + (a - b)**2*log(tan(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_204():
    f = (a + b*tan(e + f*x)**2)**2*tan(e + f*x)**6
    F = b**2*tan(e + f*x)**9/(9*f) + b*(2*a - b)*tan(e + f*x)**7/(7*f) - x*(a - b)**2 + (a - b)**2*tan(e + f*x)**5/(5*f) - (a - b)**2*tan(e + f*x)**3/(3*f) + (a - b)**2*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_205():
    f = (a + b*tan(e + f*x)**2)**2*tan(e + f*x)**4
    F = b**2*tan(e + f*x)**7/(7*f) + b*(2*a - b)*tan(e + f*x)**5/(5*f) + x*(a - b)**2 + (a - b)**2*tan(e + f*x)**3/(3*f) - (a - b)**2*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_206():
    f = (a + b*tan(e + f*x)**2)**2*tan(e + f*x)**2
    F = b**2*tan(e + f*x)**5/(5*f) + b*(2*a - b)*tan(e + f*x)**3/(3*f) - x*(a - b)**2 + (a - b)**2*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_207():
    f = (a + b*tan(e + f*x)**2)**2
    F = b**2*tan(e + f*x)**3/(3*f) + b*(2*a - b)*tan(e + f*x)/f + x*(a - b)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_208():
    f = (a + b*tan(e + f*x)**2)**2*cot(e + f*x)**2
    F = -a**2*cot(e + f*x)/f + b**2*tan(e + f*x)/f - x*(a - b)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_209():
    f = (a + b*tan(e + f*x)**2)**2*cot(e + f*x)**4
    F = -a**2*cot(e + f*x)**3/(3*f) + a*(a - 2*b)*cot(e + f*x)/f + x*(a - b)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_210():
    f = (a + b*tan(e + f*x)**2)**2*cot(e + f*x)**6
    F = -a**2*cot(e + f*x)**5/(5*f) + a*(a - 2*b)*cot(e + f*x)**3/(3*f) - x*(a - b)**2 - (a - b)**2*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_211():
    f = tan(e + f*x)**5/(a + b*tan(e + f*x)**2)
    F = -a**2*log(a + b*tan(e + f*x)**2)/(b**2*f*(2*a - 2*b)) - log(cos(e + f*x))/(f*(a - b)) + tan(e + f*x)**2/(2*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_212():
    f = tan(e + f*x)**3/(a + b*tan(e + f*x)**2)
    F = a*log(a + b*tan(e + f*x)**2)/(b*f*(2*a - 2*b)) + log(cos(e + f*x))/(f*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_213():
    f = tan(e + f*x)/(a + b*tan(e + f*x)**2)
    F = -log(a*cos(e + f*x)**2 + b*sin(e + f*x)**2)/(f*(2*a - 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_214():
    f = cot(e + f*x)/(a + b*tan(e + f*x)**2)
    F = log(cos(e + f*x))/(f*(a - b)) + b*log(a + b*tan(e + f*x)**2)/(2*a*f*(a - b)) + log(tan(e + f*x))/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_215():
    f = cot(e + f*x)**3/(a + b*tan(e + f*x)**2)
    F = -log(cos(e + f*x))/(f*(a - b)) - cot(e + f*x)**2/(2*a*f) - b**2*log(a + b*tan(e + f*x)**2)/(2*a**2*f*(a - b)) - (a + b)*log(tan(e + f*x))/(a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_216():
    f = cot(e + f*x)**5/(a + b*tan(e + f*x)**2)
    F = log(cos(e + f*x))/(f*(a - b)) - cot(e + f*x)**4/(4*a*f) + (a + b)*cot(e + f*x)**2/(2*a**2*f) + b**3*log(a + b*tan(e + f*x)**2)/(2*a**3*f*(a - b)) + (a**2 + a*b + b**2)*log(tan(e + f*x))/(a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_217():
    f = tan(e + f*x)**6/(a + b*tan(e + f*x)**2)
    F = a**(sympy.S(5)/2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(b**(sympy.S(5)/2)*f*(a - b)) - x/(a - b) + tan(e + f*x)**3/(3*b*f) - (a + b)*tan(e + f*x)/(b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_218():
    f = tan(e + f*x)**4/(a + b*tan(e + f*x)**2)
    F = -a**(sympy.S(3)/2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(b**(sympy.S(3)/2)*f*(a - b)) + x/(a - b) + tan(e + f*x)/(b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_219():
    f = tan(e + f*x)**2/(a + b*tan(e + f*x)**2)
    F = sqrt(a)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(sqrt(b)*f*(a - b)) - x/(a - b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_220():
    f = 1/(a + b*tan(e + f*x)**2)
    F = x/(a - b) - sqrt(b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(sqrt(a)*f*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_221():
    f = cot(e + f*x)**2/(a + b*tan(e + f*x)**2)
    F = -x/(a - b) - cot(e + f*x)/(a*f) + b**(sympy.S(3)/2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(a**(sympy.S(3)/2)*f*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_222():
    f = cot(e + f*x)**4/(a + b*tan(e + f*x)**2)
    F = x/(a - b) - cot(e + f*x)**3/(3*a*f) + (a + b)*cot(e + f*x)/(a**2*f) - b**(sympy.S(5)/2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(a**(sympy.S(5)/2)*f*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_223():
    f = cot(e + f*x)**6/(a + b*tan(e + f*x)**2)
    F = -x/(a - b) - cot(e + f*x)**5/(5*a*f) + (a + b)*cot(e + f*x)**3/(3*a**2*f) - (a**2 + a*b + b**2)*cot(e + f*x)/(a**3*f) + b**(sympy.S(7)/2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(a**(sympy.S(7)/2)*f*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_224():
    f = tan(e + f*x)**5/(a + b*tan(e + f*x)**2)**2
    F = a**2/(b**2*f*(a + b*tan(e + f*x)**2)*(2*a - 2*b)) + a*(a - 2*b)*log(a + b*tan(e + f*x)**2)/(2*b**2*f*(a - b)**2) - log(cos(e + f*x))/(f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_225():
    f = tan(e + f*x)**3/(a + b*tan(e + f*x)**2)**2
    F = -a/(b*f*(a + b*tan(e + f*x)**2)*(2*a - 2*b)) + log(a*cos(e + f*x)**2 + b*sin(e + f*x)**2)/(2*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_226():
    f = tan(e + f*x)/(a + b*tan(e + f*x)**2)**2
    F = 1/(f*(a + b*tan(e + f*x)**2)*(2*a - 2*b)) - log(a*cos(e + f*x)**2 + b*sin(e + f*x)**2)/(2*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_227():
    f = cot(e + f*x)/(a + b*tan(e + f*x)**2)**2
    F = log(cos(e + f*x))/(f*(a - b)**2) - b/(2*a*f*(a - b)*(a + b*tan(e + f*x)**2)) + b*(2*a - b)*log(a + b*tan(e + f*x)**2)/(2*a**2*f*(a - b)**2) + log(tan(e + f*x))/(a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_228():
    f = cot(e + f*x)**3/(a + b*tan(e + f*x)**2)**2
    F = -log(cos(e + f*x))/(f*(a - b)**2) + b**2/(2*a**2*f*(a - b)*(a + b*tan(e + f*x)**2)) - cot(e + f*x)**2/(2*a**2*f) - b**2*(3*a - 2*b)*log(a + b*tan(e + f*x)**2)/(2*a**3*f*(a - b)**2) - (a + 2*b)*log(tan(e + f*x))/(a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_229():
    f = cot(e + f*x)**5/(a + b*tan(e + f*x)**2)**2
    F = log(cos(e + f*x))/(f*(a - b)**2) - cot(e + f*x)**4/(4*a**2*f) - b**3/(2*a**3*f*(a - b)*(a + b*tan(e + f*x)**2)) + (a + 2*b)*cot(e + f*x)**2/(2*a**3*f) + b**3*(4*a - 3*b)*log(a + b*tan(e + f*x)**2)/(2*a**4*f*(a - b)**2) + (a**2 + 2*a*b + 3*b**2)*log(tan(e + f*x))/(a**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_230():
    f = tan(e + f*x)**6/(a + b*tan(e + f*x)**2)**2
    F = -a**(sympy.S(3)/2)*(3*a - 5*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(2*b**(sympy.S(5)/2)*f*(a - b)**2) - a*tan(e + f*x)**3/(b*f*(a + b*tan(e + f*x)**2)*(2*a - 2*b)) - x/(a - b)**2 + (3*a - 2*b)*tan(e + f*x)/(b**2*f*(2*a - 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_231():
    f = tan(e + f*x)**4/(a + b*tan(e + f*x)**2)**2
    F = sqrt(a)*(a - 3*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(2*b**(sympy.S(3)/2)*f*(a - b)**2) - a*tan(e + f*x)/(b*f*(a + b*tan(e + f*x)**2)*(2*a - 2*b)) + x/(a - b)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_232():
    f = tan(e + f*x)**2/(a + b*tan(e + f*x)**2)**2
    F = -x/(a - b)**2 + tan(e + f*x)/(f*(a + b*tan(e + f*x)**2)*(2*a - 2*b)) + (a + b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(2*sqrt(a)*sqrt(b)*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_233():
    f = (a + b*tan(e + f*x)**2)**(-2)
    F = x/(a - b)**2 - b*tan(e + f*x)/(2*a*f*(a - b)*(a + b*tan(e + f*x)**2)) - sqrt(b)*(3*a - b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_234():
    f = cot(e + f*x)**2/(a + b*tan(e + f*x)**2)**2
    F = -x/(a - b)**2 - b*cot(e + f*x)/(2*a*f*(a - b)*(a + b*tan(e + f*x)**2)) - (2*a - 3*b)*cot(e + f*x)/(2*a**2*f*(a - b)) + b**(sympy.S(3)/2)*(5*a - 3*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(2*a**(sympy.S(5)/2)*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_235():
    f = cot(e + f*x)**4/(a + b*tan(e + f*x)**2)**2
    F = x/(a - b)**2 - b*cot(e + f*x)**3/(2*a*f*(a - b)*(a + b*tan(e + f*x)**2)) - (2*a - 5*b)*cot(e + f*x)**3/(6*a**2*f*(a - b)) + (2*a**2 + 2*a*b - 5*b**2)*cot(e + f*x)/(2*a**3*f*(a - b)) - b**(sympy.S(5)/2)*(7*a - 5*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(2*a**(sympy.S(7)/2)*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_236():
    f = cot(e + f*x)**6/(a + b*tan(e + f*x)**2)**2
    F = -x/(a - b)**2 - b*cot(e + f*x)**5/(2*a*f*(a - b)*(a + b*tan(e + f*x)**2)) - (2*a - 7*b)*cot(e + f*x)**5/(10*a**2*f*(a - b)) + (2*a**2 + 2*a*b - 7*b**2)*cot(e + f*x)**3/(6*a**3*f*(a - b)) - (2*a**3 + 2*a**2*b + 2*a*b**2 - 7*b**3)*cot(e + f*x)/(2*a**4*f*(a - b)) + b**(sympy.S(7)/2)*(9*a - 7*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(2*a**(sympy.S(9)/2)*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_237():
    f = tan(e + f*x)**5/(a + b*tan(e + f*x)**2)**3
    F = a**2/(b**2*f*(a + b*tan(e + f*x)**2)**2*(4*a - 4*b)) - a*(a - 2*b)/(2*b**2*f*(a - b)**2*(a + b*tan(e + f*x)**2)) - log(a*cos(e + f*x)**2 + b*sin(e + f*x)**2)/(2*f*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_238():
    f = tan(e + f*x)**3/(a + b*tan(e + f*x)**2)**3
    F = -a/(b*f*(a + b*tan(e + f*x)**2)**2*(4*a - 4*b)) - 1/(2*f*(a - b)**2*(a + b*tan(e + f*x)**2)) + log(a*cos(e + f*x)**2 + b*sin(e + f*x)**2)/(2*f*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_239():
    f = tan(e + f*x)/(a + b*tan(e + f*x)**2)**3
    F = 1/(f*(a + b*tan(e + f*x)**2)**2*(4*a - 4*b)) + 1/(2*f*(a - b)**2*(a + b*tan(e + f*x)**2)) - log(a*cos(e + f*x)**2 + b*sin(e + f*x)**2)/(2*f*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_240():
    f = cot(e + f*x)/(a + b*tan(e + f*x)**2)**3
    F = log(cos(e + f*x))/(f*(a - b)**3) - b/(4*a*f*(a - b)*(a + b*tan(e + f*x)**2)**2) - b*(2*a - b)/(2*a**2*f*(a - b)**2*(a + b*tan(e + f*x)**2)) + b*(3*a**2 - 3*a*b + b**2)*log(a + b*tan(e + f*x)**2)/(2*a**3*f*(a - b)**3) + log(tan(e + f*x))/(a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_241():
    f = cot(e + f*x)**3/(a + b*tan(e + f*x)**2)**3
    F = -log(cos(e + f*x))/(f*(a - b)**3) + b**2/(4*a**2*f*(a - b)*(a + b*tan(e + f*x)**2)**2) + b**2*(3*a - 2*b)/(2*a**3*f*(a - b)**2*(a + b*tan(e + f*x)**2)) - cot(e + f*x)**2/(2*a**3*f) - b**2*(6*a**2 - 8*a*b + 3*b**2)*log(a + b*tan(e + f*x)**2)/(2*a**4*f*(a - b)**3) - (a + 3*b)*log(tan(e + f*x))/(a**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_242():
    f = cot(e + f*x)**5/(a + b*tan(e + f*x)**2)**3
    F = log(cos(e + f*x))/(f*(a - b)**3) - b**3/(4*a**3*f*(a - b)*(a + b*tan(e + f*x)**2)**2) - cot(e + f*x)**4/(4*a**3*f) - b**3*(4*a - 3*b)/(2*a**4*f*(a - b)**2*(a + b*tan(e + f*x)**2)) + (a + 3*b)*cot(e + f*x)**2/(2*a**4*f) + b**3*(10*a**2 - 15*a*b + 6*b**2)*log(a + b*tan(e + f*x)**2)/(2*a**5*f*(a - b)**3) + (a**2 + 3*a*b + 6*b**2)*log(tan(e + f*x))/(a**5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_243():
    f = tan(e + f*x)**6/(a + b*tan(e + f*x)**2)**3
    F = sqrt(a)*(3*a**2 - 10*a*b + 15*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(8*b**(sympy.S(5)/2)*f*(a - b)**3) - a*tan(e + f*x)**3/(b*f*(a + b*tan(e + f*x)**2)**2*(4*a - 4*b)) - a*(3*a - 7*b)*tan(e + f*x)/(8*b**2*f*(a - b)**2*(a + b*tan(e + f*x)**2)) - x/(a - b)**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_244():
    f = tan(e + f*x)**4/(a + b*tan(e + f*x)**2)**3
    F = -a*tan(e + f*x)/(b*f*(a + b*tan(e + f*x)**2)**2*(4*a - 4*b)) + x/(a - b)**3 + (a - 5*b)*tan(e + f*x)/(8*b*f*(a - b)**2*(a + b*tan(e + f*x)**2)) + (a**2 - 6*a*b - 3*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(8*sqrt(a)*b**(sympy.S(3)/2)*f*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_245():
    f = tan(e + f*x)**2/(a + b*tan(e + f*x)**2)**3
    F = -x/(a - b)**3 + tan(e + f*x)/(f*(a + b*tan(e + f*x)**2)**2*(4*a - 4*b)) + (3*a + b)*tan(e + f*x)/(8*a*f*(a - b)**2*(a + b*tan(e + f*x)**2)) + (3*a**2 + 6*a*b - b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(8*a**(sympy.S(3)/2)*sqrt(b)*f*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_246():
    f = (a + b*tan(e + f*x)**2)**(-3)
    F = x/(a - b)**3 - b*tan(e + f*x)/(4*a*f*(a - b)*(a + b*tan(e + f*x)**2)**2) - b*(7*a - 3*b)*tan(e + f*x)/(8*a**2*f*(a - b)**2*(a + b*tan(e + f*x)**2)) - sqrt(b)*(15*a**2 - 10*a*b + 3*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*f*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_247():
    f = cot(e + f*x)**2/(a + b*tan(e + f*x)**2)**3
    F = -x/(a - b)**3 - b*cot(e + f*x)/(4*a*f*(a - b)*(a + b*tan(e + f*x)**2)**2) - b*(9*a - 5*b)*cot(e + f*x)/(8*a**2*f*(a - b)**2*(a + b*tan(e + f*x)**2)) - (8*a**2 - 27*a*b + 15*b**2)*cot(e + f*x)/(8*a**3*f*(a - b)**2) + b**(sympy.S(3)/2)*(35*a**2 - 42*a*b + 15*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(8*a**(sympy.S(7)/2)*f*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_248():
    f = cot(e + f*x)**4/(a + b*tan(e + f*x)**2)**3
    F = x/(a - b)**3 - b*cot(e + f*x)**3/(4*a*f*(a - b)*(a + b*tan(e + f*x)**2)**2) - b*(11*a - 7*b)*cot(e + f*x)**3/(8*a**2*f*(a - b)**2*(a + b*tan(e + f*x)**2)) - (8*a**2 - 55*a*b + 35*b**2)*cot(e + f*x)**3/(24*a**3*f*(a - b)**2) + (8*a**3 + 8*a**2*b - 55*a*b**2 + 35*b**3)*cot(e + f*x)/(8*a**4*f*(a - b)**2) - b**(sympy.S(5)/2)*(63*a**2 - 90*a*b + 35*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(8*a**(sympy.S(9)/2)*f*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_249():
    f = cot(e + f*x)**6/(a + b*tan(e + f*x)**2)**3
    F = -x/(a - b)**3 - b*cot(e + f*x)**5/(4*a*f*(a - b)*(a + b*tan(e + f*x)**2)**2) - b*(13*a - 9*b)*cot(e + f*x)**5/(8*a**2*f*(a - b)**2*(a + b*tan(e + f*x)**2)) - (8*a**2 - 91*a*b + 63*b**2)*cot(e + f*x)**5/(40*a**3*f*(a - b)**2) + (8*a**3 + 8*a**2*b - 91*a*b**2 + 63*b**3)*cot(e + f*x)**3/(24*a**4*f*(a - b)**2) - (8*a**4 + 8*a**3*b + 8*a**2*b**2 - 91*a*b**3 + 63*b**4)*cot(e + f*x)/(8*a**5*f*(a - b)**2) + b**(sympy.S(7)/2)*(99*a**2 - 154*a*b + 63*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a))/(8*a**(sympy.S(11)/2)*f*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_250():
    f = (a + b*tan(c + d*x)**2)**4
    F = b**4*tan(c + d*x)**7/(7*d) + b**3*(4*a - b)*tan(c + d*x)**5/(5*d) + b**2*(6*a**2 - 4*a*b + b**2)*tan(c + d*x)**3/(3*d) + b*(2*a - b)*(2*a**2 - 2*a*b + b**2)*tan(c + d*x)/d + x*(a - b)**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_251():
    f = (a + b*tan(c + d*x)**2)**3
    F = b**3*tan(c + d*x)**5/(5*d) + b**2*(3*a - b)*tan(c + d*x)**3/(3*d) + b*(3*a**2 - 3*a*b + b**2)*tan(c + d*x)/d + x*(a - b)**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_252():
    f = (a + b*tan(c + d*x)**2)**2
    F = b**2*tan(c + d*x)**3/(3*d) + b*(2*a - b)*tan(c + d*x)/d + x*(a - b)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_253():
    f = a + b*tan(c + d*x)**2
    F = a*x - b*x + b*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_254():
    f = 1/(a + b*tan(c + d*x)**2)
    F = x/(a - b) - sqrt(b)*atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(sqrt(a)*d*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_255():
    f = (a + b*tan(c + d*x)**2)**(-2)
    F = x/(a - b)**2 - b*tan(c + d*x)/(2*a*d*(a - b)*(a + b*tan(c + d*x)**2)) - sqrt(b)*(3*a - b)*atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_256():
    f = (a + b*tan(c + d*x)**2)**(-3)
    F = x/(a - b)**3 - b*tan(c + d*x)/(4*a*d*(a - b)*(a + b*tan(c + d*x)**2)**2) - b*(7*a - 3*b)*tan(c + d*x)/(8*a**2*d*(a - b)**2*(a + b*tan(c + d*x)**2)) - sqrt(b)*(15*a**2 - 10*a*b + 3*b**2)*atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_257():
    f = sqrt(a*tan(x)**2 + a)*tan(x)**4
    F = 3*sqrt(a*sec(x)**2)*cos(x)*atanh(sin(x))/8 + sqrt(a*sec(x)**2)*tan(x)**3/4 - 3*sqrt(a*sec(x)**2)*tan(x)/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_258():
    f = sqrt(a*tan(x)**2 + a)*tan(x)**3
    F = -sqrt(a*sec(x)**2) + (a*sec(x)**2)**(sympy.S(3)/2)/(3*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_259():
    f = sqrt(a*tan(x)**2 + a)*tan(x)**2
    F = -sqrt(a*sec(x)**2)*cos(x)*atanh(sin(x))/2 + sqrt(a*sec(x)**2)*tan(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_260():
    f = sqrt(a*tan(x)**2 + a)*tan(x)
    F = sqrt(a*sec(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_261():
    f = sqrt(a*tan(x)**2 + a)*cot(x)
    F = -sqrt(a)*atanh(sqrt(a*sec(x)**2)/sqrt(a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_262():
    f = sqrt(a*tan(x)**2 + a)*cot(x)**2
    F = -sqrt(a*sec(x)**2)*cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_263():
    f = sqrt(a*tan(x)**2 + a)*cot(x)**3
    F = sqrt(a)*atanh(sqrt(a*sec(x)**2)/sqrt(a))/2 - sqrt(a*sec(x)**2)*cot(x)**2/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_264():
    f = sqrt(a*tan(x)**2 + a)*cot(x)**4
    F = -sqrt(a*sec(x)**2)*cot(x)*csc(x)**2/3 + sqrt(a*sec(x)**2)*cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_265():
    f = sqrt(a*tan(c + d*x)**2 + a)
    F = sqrt(a)*atanh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x)**2))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_266():
    f = (a*tan(x)**2 + a)**(sympy.S(3)/2)*tan(x)**3
    F = -(a*sec(x)**2)**(sympy.S(3)/2)/3 + (a*sec(x)**2)**(sympy.S(5)/2)/(5*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_267():
    f = (a*tan(x)**2 + a)**(sympy.S(3)/2)*tan(x)**2
    F = -a*sqrt(a*sec(x)**2)*cos(x)*atanh(sin(x))/8 + a*sqrt(a*sec(x)**2)*tan(x)*sec(x)**2/4 - a*sqrt(a*sec(x)**2)*tan(x)/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_268():
    f = (a*tan(x)**2 + a)**(sympy.S(3)/2)*tan(x)
    F = (a*sec(x)**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_269():
    f = (a*tan(x)**2 + a)**(sympy.S(3)/2)*cot(x)
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a*sec(x)**2)/sqrt(a)) + a*sqrt(a*sec(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_270():
    f = (a*tan(x)**2 + a)**(sympy.S(3)/2)*cot(x)**2
    F = a*sqrt(a*sec(x)**2)*cos(x)*atanh(sin(x)) - a*sqrt(a*sec(x)**2)*cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_271():
    f = (a*tan(c + d*x)**2 + a)**(sympy.S(3)/2)
    F = a**(sympy.S(3)/2)*atanh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x)**2))/(2*d) + a*sqrt(a*sec(c + d*x)**2)*tan(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_272():
    f = (a*tan(c + d*x)**2 + a)**(sympy.S(5)/2)
    F = 3*a**(sympy.S(5)/2)*atanh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x)**2))/(8*d) + 3*a**2*sqrt(a*sec(c + d*x)**2)*tan(c + d*x)/(8*d) + a*(a*sec(c + d*x)**2)**(sympy.S(3)/2)*tan(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_273():
    f = tan(x)**3/sqrt(a*tan(x)**2 + a)
    F = 1/sqrt(a*sec(x)**2) + sqrt(a*sec(x)**2)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_274():
    f = tan(x)**2/sqrt(a*tan(x)**2 + a)
    F = -tan(x)/sqrt(a*sec(x)**2) + atanh(sin(x))*sec(x)/sqrt(a*sec(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_275():
    f = tan(x)/sqrt(a*tan(x)**2 + a)
    F = -1/sqrt(a*sec(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_276():
    f = cot(x)/sqrt(a*tan(x)**2 + a)
    F = 1/sqrt(a*sec(x)**2) - atanh(sqrt(a*sec(x)**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_277():
    f = cot(x)**2/sqrt(a*tan(x)**2 + a)
    F = -tan(x)/sqrt(a*sec(x)**2) - csc(x)*sec(x)/sqrt(a*sec(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_278():
    f = tan(x)**3/(a*tan(x)**2 + a)**(sympy.S(3)/2)
    F = 1/(3*(a*sec(x)**2)**(sympy.S(3)/2)) - 1/(a*sqrt(a*sec(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_279():
    f = tan(x)**2/(a*tan(x)**2 + a)**(sympy.S(3)/2)
    F = sin(x)**2*tan(x)/(3*a*sqrt(a*sec(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_280():
    f = tan(x)/(a*tan(x)**2 + a)**(sympy.S(3)/2)
    F = -1/(3*(a*sec(x)**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_281():
    f = cot(x)/(a*tan(x)**2 + a)**(sympy.S(3)/2)
    F = 1/(3*(a*sec(x)**2)**(sympy.S(3)/2)) + 1/(a*sqrt(a*sec(x)**2)) - atanh(sqrt(a*sec(x)**2)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_282():
    f = cot(x)**2/(a*tan(x)**2 + a)**(sympy.S(3)/2)
    F = sin(x)**2*tan(x)/(3*a*sqrt(a*sec(x)**2)) - 2*tan(x)/(a*sqrt(a*sec(x)**2)) - csc(x)*sec(x)/(a*sqrt(a*sec(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_283():
    f = 1/sqrt(a*tan(c + d*x)**2 + a)
    F = tan(c + d*x)/(d*sqrt(a*sec(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_284():
    f = (a*tan(c + d*x)**2 + a)**(sympy.S(-3)/2)
    F = tan(c + d*x)/(3*d*(a*sec(c + d*x)**2)**(sympy.S(3)/2)) + 2*tan(c + d*x)/(3*a*d*sqrt(a*sec(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_285():
    f = (a*tan(c + d*x)**2 + a)**(sympy.S(-5)/2)
    F = tan(c + d*x)/(5*d*(a*sec(c + d*x)**2)**(sympy.S(5)/2)) + 4*tan(c + d*x)/(15*a*d*(a*sec(c + d*x)**2)**(sympy.S(3)/2)) + 8*tan(c + d*x)/(15*a**2*d*sqrt(a*sec(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_286():
    f = (a*tan(c + d*x)**2 + a)**(sympy.S(-7)/2)
    F = tan(c + d*x)/(7*d*(a*sec(c + d*x)**2)**(sympy.S(7)/2)) + 6*tan(c + d*x)/(35*a*d*(a*sec(c + d*x)**2)**(sympy.S(5)/2)) + 8*tan(c + d*x)/(35*a**2*d*(a*sec(c + d*x)**2)**(sympy.S(3)/2)) + 16*tan(c + d*x)/(35*a**3*d*sqrt(a*sec(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_287():
    f = (tan(x)**2 + 1)**(sympy.S(3)/2)
    F = sqrt(sec(x)**2)*tan(x)/2 + asinh(tan(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_288():
    f = sqrt(tan(x)**2 + 1)
    F = asinh(tan(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_289():
    f = 1/sqrt(tan(x)**2 + 1)
    F = tan(x)/sqrt(sec(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_290():
    f = (-tan(x)**2 - 1)**(sympy.S(3)/2)
    F = -sqrt(-sec(x)**2)*tan(x)/2 + atan(tan(x)/sqrt(-sec(x)**2))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_291():
    f = sqrt(-tan(x)**2 - 1)
    F = -atan(tan(x)/sqrt(-sec(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_292():
    f = 1/sqrt(-tan(x)**2 - 1)
    F = tan(x)/sqrt(-sec(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_293():
    f = sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)**5
    F = -sqrt(a - b)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/f + sqrt(a + b*tan(e + f*x)**2)/f - (a + b)*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)/(3*b**2*f) + (a + b*tan(e + f*x)**2)**(sympy.S(5)/2)/(5*b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_294():
    f = sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)**3
    F = sqrt(a - b)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/f - sqrt(a + b*tan(e + f*x)**2)/f + (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)/(3*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_295():
    f = sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)
    F = -sqrt(a - b)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/f + sqrt(a + b*tan(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_296():
    f = sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)
    F = -sqrt(a)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a))/f + sqrt(a - b)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_297():
    f = sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**3
    F = -sqrt(a - b)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/f - sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**2/(2*f) + (2*a - b)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a))/(2*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_298():
    f = sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**5
    F = sqrt(a - b)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/f - sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**4/(4*f) + sqrt(a + b*tan(e + f*x)**2)*(4*a - b)*cot(e + f*x)**2/(8*a*f) - (8*a**2 - 4*a*b - b**2)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a))/(8*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_299():
    f = sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)**6
    F = -sqrt(a - b)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f + sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)**5/(6*f) + (a - 6*b)*sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)**3/(24*b*f) - (a - 2*b)*(a + 4*b)*sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)/(16*b**2*f) + (a**3 + 2*a**2*b + 8*a*b**2 - 16*b**3)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(16*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_300():
    f = sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)**4
    F = sqrt(a - b)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f + sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)**3/(4*f) + (a - 4*b)*sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)/(8*b*f) - (a**2 + 4*a*b - 8*b**2)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(8*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_301():
    f = sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)**2
    F = -sqrt(a - b)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f + sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)/(2*f) + (a - 2*b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(2*sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_302():
    f = sqrt(a + b*tan(e + f*x)**2)
    F = sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f + sqrt(a - b)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_303():
    f = sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**2
    F = -sqrt(a - b)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f - sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_304():
    f = sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**4
    F = sqrt(a - b)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f - sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**3/(3*f) + sqrt(a + b*tan(e + f*x)**2)*(3*a - b)*cot(e + f*x)/(3*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_305():
    f = sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**6
    F = -sqrt(a - b)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f - sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**5/(5*f) + sqrt(a + b*tan(e + f*x)**2)*(5*a - b)*cot(e + f*x)**3/(15*a*f) - sqrt(a + b*tan(e + f*x)**2)*(15*a**2 - 5*a*b - 2*b**2)*cot(e + f*x)/(15*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_306():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)**5
    F = -(a - b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/f + (a - b)*sqrt(a + b*tan(e + f*x)**2)/f + (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)/(3*f) - (a + b)*(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)/(5*b**2*f) + (a + b*tan(e + f*x)**2)**(sympy.S(7)/2)/(7*b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_307():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)**3
    F = (a - b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/f - (a - b)*sqrt(a + b*tan(e + f*x)**2)/f - (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)/(3*f) + (a + b*tan(e + f*x)**2)**(sympy.S(5)/2)/(5*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_308():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)
    F = -(a - b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/f + (a - b)*sqrt(a + b*tan(e + f*x)**2)/f + (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_309():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a))/f + b*sqrt(a + b*tan(e + f*x)**2)/f + (a - b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_310():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**3
    F = sqrt(a)*(2*a - 3*b)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a))/(2*f) - a*sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**2/(2*f) - (a - b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_311():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**5
    F = -a*sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**4/(4*f) + (a - b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/f + sqrt(a + b*tan(e + f*x)**2)*(4*a - 5*b)*cot(e + f*x)**2/(8*f) - (8*a**2 - 12*a*b + 3*b**2)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a))/(8*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_312():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)**6
    F = b*sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)**7/(8*f) - (a - b)**(sympy.S(3)/2)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f + sqrt(a + b*tan(e + f*x)**2)*(9*a - 8*b)*tan(e + f*x)**5/(48*f) + sqrt(a + b*tan(e + f*x)**2)*(3*a**2 - 56*a*b + 48*b**2)*tan(e + f*x)**3/(192*b*f) - sqrt(a + b*tan(e + f*x)**2)*(3*a**3 + 8*a**2*b - 80*a*b**2 + 64*b**3)*tan(e + f*x)/(128*b**2*f) + (3*a**4 + 8*a**3*b + 48*a**2*b**2 - 192*a*b**3 + 128*b**4)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(128*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_313():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)**4
    F = b*sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)**5/(6*f) + (a - b)**(sympy.S(3)/2)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f + sqrt(a + b*tan(e + f*x)**2)*(7*a - 6*b)*tan(e + f*x)**3/(24*f) + sqrt(a + b*tan(e + f*x)**2)*(a**2 - 10*a*b + 8*b**2)*tan(e + f*x)/(16*b*f) - (a**3 + 6*a**2*b - 24*a*b**2 + 16*b**3)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(16*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_314():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)**2
    F = b*sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)**3/(4*f) - (a - b)**(sympy.S(3)/2)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f + sqrt(a + b*tan(e + f*x)**2)*(5*a - 4*b)*tan(e + f*x)/(8*f) + (3*a**2 - 12*a*b + 8*b**2)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(8*sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_315():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = sqrt(b)*(3*a - 2*b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(2*f) + b*sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)/(2*f) + (a - b)**(sympy.S(3)/2)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_316():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**2
    F = -a*sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)/f + b**(sympy.S(3)/2)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f - (a - b)**(sympy.S(3)/2)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_317():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**4
    F = -a*sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**3/(3*f) + (a - b)**(sympy.S(3)/2)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f + sqrt(a + b*tan(e + f*x)**2)*(3*a - 4*b)*cot(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_318():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**6
    F = -a*sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**5/(5*f) - (a - b)**(sympy.S(3)/2)*atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/f + sqrt(a + b*tan(e + f*x)**2)*(5*a - 6*b)*cot(e + f*x)**3/(15*f) - sqrt(a + b*tan(e + f*x)**2)*(15*a**2 - 20*a*b + 3*b**2)*cot(e + f*x)/(15*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_319():
    f = (a + b*tan(c + d*x)**2)**(sympy.S(5)/2)
    F = sqrt(b)*(15*a**2 - 20*a*b + 8*b**2)*atanh(sqrt(b)*tan(c + d*x)/sqrt(a + b*tan(c + d*x)**2))/(8*d) + b*(a + b*tan(c + d*x)**2)**(sympy.S(3)/2)*tan(c + d*x)/(4*d) + b*sqrt(a + b*tan(c + d*x)**2)*(7*a - 4*b)*tan(c + d*x)/(8*d) + (a - b)**(sympy.S(5)/2)*atan(sqrt(a - b)*tan(c + d*x)/sqrt(a + b*tan(c + d*x)**2))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_320():
    f = tan(e + f*x)**5/sqrt(a + b*tan(e + f*x)**2)
    F = -atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*sqrt(a - b)) - (a + b)*sqrt(a + b*tan(e + f*x)**2)/(b**2*f) + (a + b*tan(e + f*x)**2)**(sympy.S(3)/2)/(3*b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_321():
    f = tan(e + f*x)**3/sqrt(a + b*tan(e + f*x)**2)
    F = atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*sqrt(a - b)) + sqrt(a + b*tan(e + f*x)**2)/(b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_322():
    f = tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2)
    F = -atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_323():
    f = cot(e + f*x)/sqrt(a + b*tan(e + f*x)**2)
    F = atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*sqrt(a - b)) - atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_324():
    f = cot(e + f*x)**3/sqrt(a + b*tan(e + f*x)**2)
    F = -atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*sqrt(a - b)) - sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**2/(2*a*f) + (2*a + b)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a))/(2*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_325():
    f = cot(e + f*x)**5/sqrt(a + b*tan(e + f*x)**2)
    F = atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*sqrt(a - b)) - sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**4/(4*a*f) + sqrt(a + b*tan(e + f*x)**2)*(4*a + 3*b)*cot(e + f*x)**2/(8*a**2*f) - (8*a**2 + 4*a*b + 3*b**2)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a))/(8*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_326():
    f = tan(e + f*x)**6/sqrt(a + b*tan(e + f*x)**2)
    F = -atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*sqrt(a - b)) + sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)**3/(4*b*f) - sqrt(a + b*tan(e + f*x)**2)*(3*a + 4*b)*tan(e + f*x)/(8*b**2*f) + (3*a**2 + 4*a*b + 8*b**2)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(8*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_327():
    f = tan(e + f*x)**4/sqrt(a + b*tan(e + f*x)**2)
    F = atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*sqrt(a - b)) + sqrt(a + b*tan(e + f*x)**2)*tan(e + f*x)/(2*b*f) - (a + 2*b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(2*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_328():
    f = tan(e + f*x)**2/sqrt(a + b*tan(e + f*x)**2)
    F = -atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*sqrt(a - b)) + atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_329():
    f = 1/sqrt(a + b*tan(e + f*x)**2)
    F = atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_330():
    f = cot(e + f*x)**2/sqrt(a + b*tan(e + f*x)**2)
    F = -atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*sqrt(a - b)) - sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_331():
    f = cot(e + f*x)**4/sqrt(a + b*tan(e + f*x)**2)
    F = atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*sqrt(a - b)) - sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**3/(3*a*f) + sqrt(a + b*tan(e + f*x)**2)*(3*a + 2*b)*cot(e + f*x)/(3*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_332():
    f = cot(e + f*x)**6/sqrt(a + b*tan(e + f*x)**2)
    F = -atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*sqrt(a - b)) - sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**5/(5*a*f) + sqrt(a + b*tan(e + f*x)**2)*(5*a + 4*b)*cot(e + f*x)**3/(15*a**2*f) - sqrt(a + b*tan(e + f*x)**2)*(15*a**2 + 10*a*b + 8*b**2)*cot(e + f*x)/(15*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_333():
    f = tan(e + f*x)**5/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = a**2/(b**2*f*(a - b)*sqrt(a + b*tan(e + f*x)**2)) - atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*(a - b)**(sympy.S(3)/2)) + sqrt(a + b*tan(e + f*x)**2)/(b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_334():
    f = tan(e + f*x)**3/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -a/(b*f*(a - b)*sqrt(a + b*tan(e + f*x)**2)) + atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_335():
    f = tan(e + f*x)/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = 1/(f*(a - b)*sqrt(a + b*tan(e + f*x)**2)) - atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_336():
    f = cot(e + f*x)/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*(a - b)**(sympy.S(3)/2)) - b/(a*f*(a - b)*sqrt(a + b*tan(e + f*x)**2)) - atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_337():
    f = cot(e + f*x)**3/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*(a - b)**(sympy.S(3)/2)) - cot(e + f*x)**2/(2*a*f*sqrt(a + b*tan(e + f*x)**2)) - b*(a - 3*b)/(2*a**2*f*(a - b)*sqrt(a + b*tan(e + f*x)**2)) + (2*a + 3*b)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a))/(2*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_338():
    f = cot(e + f*x)**5/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*(a - b)**(sympy.S(3)/2)) - cot(e + f*x)**4/(4*a*f*sqrt(a + b*tan(e + f*x)**2)) + (4*a + 5*b)*cot(e + f*x)**2/(8*a**2*f*sqrt(a + b*tan(e + f*x)**2)) + b*(4*a**2 + 3*a*b - 15*b**2)/(8*a**3*f*(a - b)*sqrt(a + b*tan(e + f*x)**2)) - (8*a**2 + 12*a*b + 15*b**2)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a))/(8*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_339():
    f = tan(e + f*x)**6/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -a*tan(e + f*x)**3/(b*f*(a - b)*sqrt(a + b*tan(e + f*x)**2)) - atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(3)/2)) + sqrt(a + b*tan(e + f*x)**2)*(3*a - b)*tan(e + f*x)/(b**2*f*(2*a - 2*b)) - (3*a + 2*b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(2*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_340():
    f = tan(e + f*x)**4/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -a*tan(e + f*x)/(b*f*(a - b)*sqrt(a + b*tan(e + f*x)**2)) + atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(3)/2)) + atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_341():
    f = tan(e + f*x)**2/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = tan(e + f*x)/(f*(a - b)*sqrt(a + b*tan(e + f*x)**2)) - atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_342():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(-3)/2)
    F = atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(3)/2)) - b*tan(e + f*x)/(a*f*(a - b)*sqrt(a + b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_343():
    f = cot(e + f*x)**2/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(3)/2)) - b*cot(e + f*x)/(a*f*(a - b)*sqrt(a + b*tan(e + f*x)**2)) - (a - 2*b)*sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)/(a**2*f*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_344():
    f = cot(e + f*x)**4/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(3)/2)) - b*cot(e + f*x)**3/(a*f*(a - b)*sqrt(a + b*tan(e + f*x)**2)) - (a - 4*b)*sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**3/(3*a**2*f*(a - b)) + (a + 2*b)*sqrt(a + b*tan(e + f*x)**2)*(3*a - 4*b)*cot(e + f*x)/(3*a**3*f*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_345():
    f = cot(e + f*x)**6/(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)
    F = -atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(3)/2)) - b*cot(e + f*x)**5/(a*f*(a - b)*sqrt(a + b*tan(e + f*x)**2)) - (a - 6*b)*sqrt(a + b*tan(e + f*x)**2)*cot(e + f*x)**5/(5*a**2*f*(a - b)) + sqrt(a + b*tan(e + f*x)**2)*(5*a**2 + 4*a*b - 24*b**2)*cot(e + f*x)**3/(15*a**3*f*(a - b)) - sqrt(a + b*tan(e + f*x)**2)*(15*a**3 + 10*a**2*b + 8*a*b**2 - 48*b**3)*cot(e + f*x)/(15*a**4*f*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_346():
    f = tan(e + f*x)**5/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = a**2/(b**2*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*(3*a - 3*b)) - a*(a - 2*b)/(b**2*f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2)) - atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_347():
    f = tan(e + f*x)**3/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -a/(b*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*(3*a - 3*b)) - 1/(f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2)) + atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_348():
    f = tan(e + f*x)/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = 1/(f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*(3*a - 3*b)) + 1/(f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2)) - atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_349():
    f = cot(e + f*x)/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*(a - b)**(sympy.S(5)/2)) - b/(3*a*f*(a - b)*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - b*(2*a - b)/(a**2*f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2)) - atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_350():
    f = cot(e + f*x)**3/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*(a - b)**(sympy.S(5)/2)) - cot(e + f*x)**2/(2*a*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - b*(3*a - 5*b)/(6*a**2*f*(a - b)*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - b*(a**2 - 8*a*b + 5*b**2)/(2*a**3*f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2)) + (2*a + 5*b)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a))/(2*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_351():
    f = cot(e + f*x)**5/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a - b))/(f*(a - b)**(sympy.S(5)/2)) - cot(e + f*x)**4/(4*a*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) + (4*a + 7*b)*cot(e + f*x)**2/(8*a**2*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) + b*(12*a**2 + 15*a*b - 35*b**2)/(24*a**3*f*(a - b)*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) + b*(4*a**3 + 3*a**2*b - 50*a*b**2 + 35*b**3)/(8*a**4*f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2)) - (8*a**2 + 20*a*b + 35*b**2)*atanh(sqrt(a + b*tan(e + f*x)**2)/sqrt(a))/(8*a**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_352():
    f = tan(e + f*x)**6/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -a*tan(e + f*x)**3/(b*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*(3*a - 3*b)) - a*(a - 2*b)*tan(e + f*x)/(b**2*f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2)) - atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(5)/2)) + atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_353():
    f = tan(e + f*x)**4/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -a*tan(e + f*x)/(b*f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*(3*a - 3*b)) + atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(5)/2)) + (a - 4*b)*tan(e + f*x)/(3*b*f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_354():
    f = tan(e + f*x)**2/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = tan(e + f*x)/(f*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)*(3*a - 3*b)) - atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(5)/2)) + (2*a + b)*tan(e + f*x)/(3*a*f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_355():
    f = (a + b*tan(e + f*x)**2)**(sympy.S(-5)/2)
    F = atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(5)/2)) - b*tan(e + f*x)/(3*a*f*(a - b)*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - b*(5*a - 2*b)*tan(e + f*x)/(3*a**2*f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_356():
    f = cot(e + f*x)**2/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(5)/2)) - b*cot(e + f*x)/(3*a*f*(a - b)*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - b*(7*a - 4*b)*cot(e + f*x)/(3*a**2*f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2)) - (a - 4*b)*sqrt(a + b*tan(e + f*x)**2)*(3*a - 2*b)*cot(e + f*x)/(3*a**3*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_357():
    f = cot(e + f*x)**4/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(5)/2)) - b*cot(e + f*x)**3/(3*a*f*(a - b)*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - b*(3*a - 2*b)*cot(e + f*x)**3/(a**2*f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2)) - sqrt(a + b*tan(e + f*x)**2)*(a**2 - 12*a*b + 8*b**2)*cot(e + f*x)**3/(3*a**3*f*(a - b)**2) + (a - 2*b)*sqrt(a + b*tan(e + f*x)**2)*(3*a**2 + 8*a*b - 8*b**2)*cot(e + f*x)/(3*a**4*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_358():
    f = cot(e + f*x)**6/(a + b*tan(e + f*x)**2)**(sympy.S(5)/2)
    F = -atan(sqrt(a - b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2))/(f*(a - b)**(sympy.S(5)/2)) - b*cot(e + f*x)**5/(3*a*f*(a - b)*(a + b*tan(e + f*x)**2)**(sympy.S(3)/2)) - b*(11*a - 8*b)*cot(e + f*x)**5/(3*a**2*f*(a - b)**2*sqrt(a + b*tan(e + f*x)**2)) - sqrt(a + b*tan(e + f*x)**2)*(a**2 - 22*a*b + 16*b**2)*cot(e + f*x)**5/(5*a**3*f*(a - b)**2) + sqrt(a + b*tan(e + f*x)**2)*(5*a**3 + 4*a**2*b - 88*a*b**2 + 64*b**3)*cot(e + f*x)**3/(15*a**4*f*(a - b)**2) - sqrt(a + b*tan(e + f*x)**2)*(15*a**4 + 10*a**3*b + 8*a**2*b**2 - 176*a*b**3 + 128*b**4)*cot(e + f*x)/(15*a**5*f*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_359():
    f = (b*tan(e + f*x)**2)**p*(d*tan(e + f*x))**m
    F = (b*tan(e + f*x)**2)**p*(d*tan(e + f*x))**m*tan(e + f*x)*hyper((1, m/2 + p + sympy.S.Half), (m/2 + p + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(m + 2*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_360():
    f = (d*tan(e + f*x))**m*(a + b*tan(e + f*x)**2)**p
    F = (d*tan(e + f*x))**(m + 1)*(a + b*tan(e + f*x)**2)**p*appellf1(m/2 + sympy.S.Half, 1, -p, m/2 + sympy.S(3)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/a)/(d*f*(1 + b*tan(e + f*x)**2/a)**p*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_361():
    f = (a + b*tan(e + f*x)**2)**p*tan(e + f*x)**5
    F = -(a + b*tan(e + f*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), (a + b*tan(e + f*x)**2)/(a - b))/(f*(2*a - 2*b)*(p + 1)) - (a + b)*(a + b*tan(e + f*x)**2)**(p + 1)/(2*b**2*f*(p + 1)) + (a + b*tan(e + f*x)**2)**(p + 2)/(2*b**2*f*(p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_362():
    f = (a + b*tan(e + f*x)**2)**p*tan(e + f*x)**3
    F = (a + b*tan(e + f*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), (a + b*tan(e + f*x)**2)/(a - b))/(f*(2*a - 2*b)*(p + 1)) + (a + b*tan(e + f*x)**2)**(p + 1)/(2*b*f*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_363():
    f = (a + b*tan(e + f*x)**2)**p*tan(e + f*x)
    F = -(a + b*tan(e + f*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), (a + b*tan(e + f*x)**2)/(a - b))/(f*(2*a - 2*b)*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_364():
    f = (a + b*tan(e + f*x)**2)**p*cot(e + f*x)
    F = (a + b*tan(e + f*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), (a + b*tan(e + f*x)**2)/(a - b))/(f*(2*a - 2*b)*(p + 1)) - (a + b*tan(e + f*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*tan(e + f*x)**2/a)/(2*a*f*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_365():
    f = (a + b*tan(e + f*x)**2)**p*cot(e + f*x)**3
    F = -(a + b*tan(e + f*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), (a + b*tan(e + f*x)**2)/(a - b))/(f*(2*a - 2*b)*(p + 1)) - (a + b*tan(e + f*x)**2)**(p + 1)*cot(e + f*x)**2/(2*a*f) + (a - b*p)*(a + b*tan(e + f*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*tan(e + f*x)**2/a)/(2*a**2*f*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_366():
    f = (a + b*tan(e + f*x)**2)**p*cot(e + f*x)**5
    F = (a + b*tan(e + f*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), (a + b*tan(e + f*x)**2)/(a - b))/(f*(2*a - 2*b)*(p + 1)) - (a + b*tan(e + f*x)**2)**(p + 1)*cot(e + f*x)**4/(4*a*f) + (a + b*tan(e + f*x)**2)**(p + 1)*(2*a - b*p + b)*cot(e + f*x)**2/(4*a**2*f) - (a + b*tan(e + f*x)**2)**(p + 1)*(2*a**2 - 2*a*b*p - b**2*p*(1 - p))*hyper((1, p + 1), (p + 2,), 1 + b*tan(e + f*x)**2/a)/(4*a**3*f*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_367():
    f = (a + b*tan(e + f*x)**2)**p*tan(e + f*x)**6
    F = (a + b*tan(e + f*x)**2)**p*tan(e + f*x)**7*appellf1(sympy.S(7)/2, 1, -p, sympy.S(9)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/a)/(7*f*(1 + b*tan(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_368():
    f = (a + b*tan(e + f*x)**2)**p*tan(e + f*x)**4
    F = (a + b*tan(e + f*x)**2)**p*tan(e + f*x)**5*appellf1(sympy.S(5)/2, 1, -p, sympy.S(7)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/a)/(5*f*(1 + b*tan(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_369():
    f = (a + b*tan(e + f*x)**2)**p*tan(e + f*x)**2
    F = (a + b*tan(e + f*x)**2)**p*tan(e + f*x)**3*appellf1(sympy.S(3)/2, 1, -p, sympy.S(5)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/a)/(3*f*(1 + b*tan(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_370():
    f = (a + b*tan(e + f*x)**2)**p
    F = (a + b*tan(e + f*x)**2)**p*tan(e + f*x)*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/a)/(f*(1 + b*tan(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_371():
    f = (a + b*tan(e + f*x)**2)**p*cot(e + f*x)**2
    F = -(a + b*tan(e + f*x)**2)**p*cot(e + f*x)*appellf1(sympy.S(-1)/2, 1, -p, sympy.S.Half, -tan(e + f*x)**2, -b*tan(e + f*x)**2/a)/(f*(1 + b*tan(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_372():
    f = (a + b*tan(e + f*x)**2)**p*cot(e + f*x)**4
    F = -(a + b*tan(e + f*x)**2)**p*cot(e + f*x)**3*appellf1(sympy.S(-3)/2, 1, -p, sympy.S(-1)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/a)/(3*f*(1 + b*tan(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_373():
    f = (a + b*tan(e + f*x)**2)**p*cot(e + f*x)**6
    F = -(a + b*tan(e + f*x)**2)**p*cot(e + f*x)**5*appellf1(sympy.S(-5)/2, 1, -p, sympy.S(-3)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/a)/(5*f*(1 + b*tan(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_374():
    f = (a + b*tan(c + d*x)**3)**4
    F = a*b**3*tan(c + d*x)**8/(2*d) - 2*a*b**3*tan(c + d*x)**6/(3*d) + a*b**3*tan(c + d*x)**4/d + 4*a*b*(a**2 - b**2)*log(cos(c + d*x))/d + 2*a*b*(a**2 - b**2)*tan(c + d*x)**2/d + b**4*tan(c + d*x)**11/(11*d) - b**4*tan(c + d*x)**9/(9*d) + b**4*tan(c + d*x)**7/(7*d) + b**2*(6*a**2 - b**2)*tan(c + d*x)**5/(5*d) - b**2*(6*a**2 - b**2)*tan(c + d*x)**3/(3*d) + b**2*(6*a**2 - b**2)*tan(c + d*x)/d + x*(a**4 - 6*a**2*b**2 + b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_375():
    f = (a + b*tan(c + d*x)**3)**3
    F = 3*a*b**2*tan(c + d*x)**5/(5*d) - a*b**2*tan(c + d*x)**3/d + 3*a*b**2*tan(c + d*x)/d + a*x*(a**2 - 3*b**2) + b**3*tan(c + d*x)**8/(8*d) - b**3*tan(c + d*x)**6/(6*d) + b**3*tan(c + d*x)**4/(4*d) + b*(3*a**2 - b**2)*log(cos(c + d*x))/d + b*(3*a**2 - b**2)*tan(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_376():
    f = (a + b*tan(c + d*x)**3)**2
    F = 2*a*b*log(cos(c + d*x))/d + a*b*tan(c + d*x)**2/d + b**2*tan(c + d*x)**5/(5*d) - b**2*tan(c + d*x)**3/(3*d) + b**2*tan(c + d*x)/d + x*(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_377():
    f = a + b*tan(c + d*x)**3
    F = a*x + b*log(cos(c + d*x))/d + b*tan(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_378():
    f = 1/(a + b*tan(c + d*x)**3)
    F = a*x/(a**2 + b**2) - b*log(a*cos(c + d*x)**3 + b*sin(c + d*x)**3)/(d*(3*a**2 + 3*b**2)) + sqrt(3)*b**(sympy.S(1)/3)*(a**(sympy.S(4)/3) - b**(sympy.S(4)/3))*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*tan(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*d*(a**2 + b**2)) + b**(sympy.S(1)/3)*(a**(sympy.S(4)/3) + b**(sympy.S(4)/3))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*tan(c + d*x))/(3*a**(sympy.S(2)/3)*d*(a**2 + b**2)) - b**(sympy.S(1)/3)*(a**(sympy.S(4)/3) + b**(sympy.S(4)/3))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*tan(c + d*x) + b**(sympy.S(2)/3)*tan(c + d*x)**2)/(6*a**(sympy.S(2)/3)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_379():
    f = (a + b*tan(c + d*x)**3)**(-2)
    F = -2*a*b*log(a*cos(c + d*x)**3 + b*sin(c + d*x)**3)/(3*d*(a**2 + b**2)**2) + x*(a**2 - b**2)/(a**2 + b**2)**2 + b*(a + (-a*tan(c + d*x) + b)*tan(c + d*x))/(3*a*d*(a + b*tan(c + d*x)**3)*(a**2 + b**2)) + sqrt(3)*b**(sympy.S(1)/3)*(-2*a**(sympy.S(2)/3)*b**(sympy.S(4)/3) + a**2 - b**2)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*tan(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*d*(a**2 + b**2)**2) + b**(sympy.S(1)/3)*(2*a**(sympy.S(2)/3)*b**(sympy.S(4)/3) + a**2 - b**2)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*tan(c + d*x))/(3*a**(sympy.S(1)/3)*d*(a**2 + b**2)**2) - b**(sympy.S(1)/3)*(2*a**(sympy.S(2)/3)*b**(sympy.S(4)/3) + a**2 - b**2)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*tan(c + d*x) + b**(sympy.S(2)/3)*tan(c + d*x)**2)/(6*a**(sympy.S(1)/3)*d*(a**2 + b**2)**2) + sqrt(3)*b**(sympy.S(1)/3)*(a**(sympy.S(4)/3) - 2*b**(sympy.S(4)/3))*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*tan(c + d*x))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*d*(a**2 + b**2)) + b**(sympy.S(1)/3)*(a**(sympy.S(4)/3) + 2*b**(sympy.S(4)/3))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*tan(c + d*x))/(9*a**(sympy.S(5)/3)*d*(a**2 + b**2)) - b**(sympy.S(1)/3)*(a**(sympy.S(4)/3) + 2*b**(sympy.S(4)/3))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*tan(c + d*x) + b**(sympy.S(2)/3)*tan(c + d*x)**2)/(18*a**(sympy.S(5)/3)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_380():
    f = 1/(tan(x)**3 + 1)
    F = x/2 + log(tan(x) + 1)/6 - log(tan(x)**2 - tan(x) + 1)/3 - log(cos(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_381():
    f = (a + b*tan(c + d*x)**4)**4
    F = b**4*tan(c + d*x)**15/(15*d) - b**4*tan(c + d*x)**13/(13*d) + b**3*(4*a + b)*tan(c + d*x)**11/(11*d) - b**3*(4*a + b)*tan(c + d*x)**9/(9*d) + b**2*(6*a**2 + 4*a*b + b**2)*tan(c + d*x)**7/(7*d) - b**2*(6*a**2 + 4*a*b + b**2)*tan(c + d*x)**5/(5*d) + b*(2*a + b)*(2*a**2 + 2*a*b + b**2)*tan(c + d*x)**3/(3*d) - b*(2*a + b)*(2*a**2 + 2*a*b + b**2)*tan(c + d*x)/d + x*(a + b)**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_382():
    f = (a + b*tan(c + d*x)**4)**3
    F = b**3*tan(c + d*x)**11/(11*d) - b**3*tan(c + d*x)**9/(9*d) + b**2*(3*a + b)*tan(c + d*x)**7/(7*d) - b**2*(3*a + b)*tan(c + d*x)**5/(5*d) + b*(3*a**2 + 3*a*b + b**2)*tan(c + d*x)**3/(3*d) - b*(3*a**2 + 3*a*b + b**2)*tan(c + d*x)/d + x*(a + b)**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_383():
    f = (a + b*tan(c + d*x)**4)**2
    F = b**2*tan(c + d*x)**7/(7*d) - b**2*tan(c + d*x)**5/(5*d) + b*(2*a + b)*tan(c + d*x)**3/(3*d) - b*(2*a + b)*tan(c + d*x)/d + x*(a + b)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_384():
    f = a + b*tan(c + d*x)**4
    F = a*x + b*x + b*tan(c + d*x)**3/(3*d) - b*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_385():
    f = 1/(a + b*tan(c + d*x)**4)
    F = x/(a + b) + sqrt(2)*b**(sympy.S(1)/4)*(sqrt(a) - sqrt(b))*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*tan(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*d*(a + b)) - sqrt(2)*b**(sympy.S(1)/4)*(sqrt(a) - sqrt(b))*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*tan(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*d*(a + b)) - sqrt(2)*b**(sympy.S(1)/4)*(sqrt(a) + sqrt(b))*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*tan(c + d*x) + sqrt(a) + sqrt(b)*tan(c + d*x)**2)/(8*a**(sympy.S(3)/4)*d*(a + b)) + sqrt(2)*b**(sympy.S(1)/4)*(sqrt(a) + sqrt(b))*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*tan(c + d*x) + sqrt(a) + sqrt(b)*tan(c + d*x)**2)/(8*a**(sympy.S(3)/4)*d*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_386():
    f = (a + b*tan(c + d*x)**4)**(-2)
    F = x/(a + b)**2 + b*(1 - tan(c + d*x)**2)*tan(c + d*x)/(4*a*d*(a + b)*(a + b*tan(c + d*x)**4)) + sqrt(2)*b**(sympy.S(1)/4)*(sqrt(a) - sqrt(b))*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*tan(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*d*(a + b)**2) - sqrt(2)*b**(sympy.S(1)/4)*(sqrt(a) - sqrt(b))*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*tan(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*d*(a + b)**2) - sqrt(2)*b**(sympy.S(1)/4)*(sqrt(a) + sqrt(b))*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*tan(c + d*x) + sqrt(a) + sqrt(b)*tan(c + d*x)**2)/(8*a**(sympy.S(3)/4)*d*(a + b)**2) + sqrt(2)*b**(sympy.S(1)/4)*(sqrt(a) + sqrt(b))*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*tan(c + d*x) + sqrt(a) + sqrt(b)*tan(c + d*x)**2)/(8*a**(sympy.S(3)/4)*d*(a + b)**2) + sqrt(2)*b**(sympy.S(1)/4)*(sqrt(a) - 3*sqrt(b))*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*tan(c + d*x)/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*d*(a + b)) - sqrt(2)*b**(sympy.S(1)/4)*(sqrt(a) - 3*sqrt(b))*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*tan(c + d*x)/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*d*(a + b)) - sqrt(2)*b**(sympy.S(1)/4)*(sqrt(a) + 3*sqrt(b))*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*tan(c + d*x) + sqrt(a) + sqrt(b)*tan(c + d*x)**2)/(32*a**(sympy.S(7)/4)*d*(a + b)) + sqrt(2)*b**(sympy.S(1)/4)*(sqrt(a) + 3*sqrt(b))*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*tan(c + d*x) + sqrt(a) + sqrt(b)*tan(c + d*x)**2)/(32*a**(sympy.S(7)/4)*d*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_387():
    f = sqrt(a + b*tan(c + d*x)**4)
    F = ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.atan(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4))))))**(Integer(-1))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt(Symbol('b')) * sympy.tan((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))) * ((Symbol('d') * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1)))) + (((sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('a') + Symbol('b')) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1)))) + (((sympy.sqrt(Symbol('a')) + sympy.sqrt(Symbol('b'))) * (Symbol('a') + Symbol('b')) * sympy.Function('EllipticPi')((Integer(-1) * (((sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_388():
    f = 1/sqrt(a + b*tan(c + d*x)**4)
    F = (sympy.atan(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4))))))**(Integer(-1)))) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1)))) + (((sympy.sqrt(Symbol('a')) + sympy.sqrt(Symbol('b'))) * sympy.Function('EllipticPi')((Integer(-1) * (((sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_389():
    f = sqrt(a + b*tan(x)**4)*tan(x)**3
    F = -(sympy.S.Half - tan(x)**2/4)*sqrt(a + b*tan(x)**4) + sqrt(a + b)*atanh((a - b*tan(x)**2)/(sqrt(a + b)*sqrt(a + b*tan(x)**4)))/2 + (a + 2*b)*atanh(sqrt(b)*tan(x)**2/sqrt(a + b*tan(x)**4))/(4*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_390():
    f = sqrt(a + b*tan(x)**4)*tan(x)
    F = -sqrt(b)*atanh(sqrt(b)*tan(x)**2/sqrt(a + b*tan(x)**4))/2 - sqrt(a + b)*atanh((a - b*tan(x)**2)/(sqrt(a + b)*sqrt(a + b*tan(x)**4)))/2 + sqrt(a + b*tan(x)**4)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_391():
    f = sqrt(a + b*tan(x)**4)*cot(x)
    F = -sqrt(a)*atanh(sqrt(a + b*tan(x)**4)/sqrt(a))/2 + sqrt(b)*atanh(sqrt(b)*tan(x)**2/sqrt(a + b*tan(x)**4))/2 + sqrt(a + b)*atanh((a - b*tan(x)**2)/(sqrt(a + b)*sqrt(a + b*tan(x)**4)))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_392():
    f = sqrt(a + b*tan(x)**4)*tan(x)**2
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.atan(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.tan(x)) * (sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4))))))**(Integer(-1))))) + ((Integer(3))**(Integer(-1)) * sympy.tan(x) * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4)))))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.tan(x) * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4)))))) * ((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan(x))**(Integer(2)))))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.tan(x)) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan(x))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan(x))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.tan(x)) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan(x))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan(x))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.tan(x)) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan(x))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan(x))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4)))))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('a') + Symbol('b')) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.tan(x)) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan(x))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan(x))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((sympy.sqrt(Symbol('a')) + sympy.sqrt(Symbol('b'))) * (Symbol('a') + Symbol('b')) * sympy.Function('EllipticPi')((Integer(-1) * (((sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.tan(x)) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan(x))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan(x))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_393():
    f = (a + b*tan(x)**4)**(sympy.S(3)/2)*tan(x)**3
    F = -(sympy.S(1)/6 - tan(x)**2/8)*(a + b*tan(x)**4)**(sympy.S(3)/2) + (a + b)**(sympy.S(3)/2)*atanh((a - b*tan(x)**2)/(sqrt(a + b)*sqrt(a + b*tan(x)**4)))/2 - sqrt(a + b*tan(x)**4)*(a/2 + b/2 - (3*a + 4*b)*tan(x)**2/16) + (3*a**2 + 12*a*b + 8*b**2)*atanh(sqrt(b)*tan(x)**2/sqrt(a + b*tan(x)**4))/(16*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_394():
    f = (a + b*tan(x)**4)**(sympy.S(3)/2)*tan(x)
    F = -sqrt(b)*(3*a + 2*b)*atanh(sqrt(b)*tan(x)**2/sqrt(a + b*tan(x)**4))/4 - (a + b)**(sympy.S(3)/2)*atanh((a - b*tan(x)**2)/(sqrt(a + b)*sqrt(a + b*tan(x)**4)))/2 + (a + b*tan(x)**4)**(sympy.S(3)/2)/6 + sqrt(a + b*tan(x)**4)*(a/2 - b*tan(x)**2/4 + b/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_395():
    f = (a + b*tan(x)**4)**(sympy.S(3)/2)*cot(x)
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(x)**4)/sqrt(a))/2 + a*sqrt(a + b*tan(x)**4)/2 + sqrt(b)*(3*a + 2*b)*atanh(sqrt(b)*tan(x)**2/sqrt(a + b*tan(x)**4))/4 + (a + b)**(sympy.S(3)/2)*atanh((a - b*tan(x)**2)/(sqrt(a + b)*sqrt(a + b*tan(x)**4)))/2 - sqrt(a + b*tan(x)**4)*(a/2 - b*tan(x)**2/4 + b/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_396():
    f = tan(x)**3/sqrt(a + b*tan(x)**4)
    F = atanh((a - b*tan(x)**2)/(sqrt(a + b)*sqrt(a + b*tan(x)**4)))/(2*sqrt(a + b)) + atanh(sqrt(b)*tan(x)**2/sqrt(a + b*tan(x)**4))/(2*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_397():
    f = tan(x)/sqrt(a + b*tan(x)**4)
    F = -atanh((a - b*tan(x)**2)/(sqrt(a + b)*sqrt(a + b*tan(x)**4)))/(2*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_398():
    f = cot(x)/sqrt(a + b*tan(x)**4)
    F = atanh((a - b*tan(x)**2)/(sqrt(a + b)*sqrt(a + b*tan(x)**4)))/(2*sqrt(a + b)) - atanh(sqrt(a + b*tan(x)**4)/sqrt(a))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_399():
    f = tan(x)**2/sqrt(a + b*tan(x)**4)
    F = (Integer(-1) * (sympy.atan(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.tan(x)) * (sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4))))))**(Integer(-1)))) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.tan(x)) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan(x))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan(x))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * (sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((sympy.sqrt(Symbol('a')) + sympy.sqrt(Symbol('b'))) * sympy.Function('EllipticPi')((Integer(-1) * (((sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.tan(x)) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan(x))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (sympy.tan(x))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.tan(x))**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_400():
    f = tan(x)**3/(a + b*tan(x)**4)**(sympy.S(3)/2)
    F = -(1 - tan(x)**2)/(sqrt(a + b*tan(x)**4)*(2*a + 2*b)) + atanh((a - b*tan(x)**2)/(sqrt(a + b)*sqrt(a + b*tan(x)**4)))/(2*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_401():
    f = tan(x)/(a + b*tan(x)**4)**(sympy.S(3)/2)
    F = -atanh((a - b*tan(x)**2)/(sqrt(a + b)*sqrt(a + b*tan(x)**4)))/(2*(a + b)**(sympy.S(3)/2)) + (a + b*tan(x)**2)/(2*a*(a + b)*sqrt(a + b*tan(x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_402():
    f = cot(x)/(a + b*tan(x)**4)**(sympy.S(3)/2)
    F = atanh((a - b*tan(x)**2)/(sqrt(a + b)*sqrt(a + b*tan(x)**4)))/(2*(a + b)**(sympy.S(3)/2)) + 1/(2*a*sqrt(a + b*tan(x)**4)) - (a + b*tan(x)**2)/(2*a*(a + b)*sqrt(a + b*tan(x)**4)) - atanh(sqrt(a + b*tan(x)**4)/sqrt(a))/(2*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_403():
    f = tan(x)**3/(a + b*tan(x)**4)**(sympy.S(5)/2)
    F = -(1 - tan(x)**2)/((a + b*tan(x)**4)**(sympy.S(3)/2)*(6*a + 6*b)) + atanh((a - b*tan(x)**2)/(sqrt(a + b)*sqrt(a + b*tan(x)**4)))/(2*(a + b)**(sympy.S(5)/2)) - (3*a + (-2*a + b)*tan(x)**2)/(6*a*(a + b)**2*sqrt(a + b*tan(x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_404():
    f = tan(x)/(a + b*tan(x)**4)**(sympy.S(5)/2)
    F = -atanh((a - b*tan(x)**2)/(sqrt(a + b)*sqrt(a + b*tan(x)**4)))/(2*(a + b)**(sympy.S(5)/2)) + (a + b*tan(x)**2)/(6*a*(a + b)*(a + b*tan(x)**4)**(sympy.S(3)/2)) + (3*a**2 + b*(5*a + 2*b)*tan(x)**2)/(6*a**2*(a + b)**2*sqrt(a + b*tan(x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_405():
    f = cot(x)/(a + b*tan(x)**4)**(sympy.S(5)/2)
    F = atanh((a - b*tan(x)**2)/(sqrt(a + b)*sqrt(a + b*tan(x)**4)))/(2*(a + b)**(sympy.S(5)/2)) + 1/(6*a*(a + b*tan(x)**4)**(sympy.S(3)/2)) - (a + b*tan(x)**2)/(6*a*(a + b)*(a + b*tan(x)**4)**(sympy.S(3)/2)) + 1/(2*a**2*sqrt(a + b*tan(x)**4)) - (3*a**2 + b*(5*a + 2*b)*tan(x)**2)/(6*a**2*(a + b)**2*sqrt(a + b*tan(x)**4)) - atanh(sqrt(a + b*tan(x)**4)/sqrt(a))/(2*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_406():
    f = (d*tan(e + f*x))**m*(a + b*sqrt(c*tan(e + f*x)))**2
    F = 4*a*b*(c*tan(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**m*hyper((1, m/2 + sympy.S(3)/4), (m/2 + sympy.S(7)/4,), -tan(e + f*x)**2)/(c*f*(2*m + 3)) + (d*tan(e + f*x))**m*(a**2 - b**2*sqrt(-c**2))*tan(e + f*x)*hyper((1, m + 1), (m + 2,), -c*tan(e + f*x)/sqrt(-c**2))/(2*f*(m + 1)) + (d*tan(e + f*x))**m*(a**2 + b**2*sqrt(-c**2))*tan(e + f*x)*hyper((1, m + 1), (m + 2,), c*tan(e + f*x)/sqrt(-c**2))/(2*f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_407():
    f = (d*tan(e + f*x))**m*(a + b*sqrt(c*tan(e + f*x)))
    F = a*(d*tan(e + f*x))**m*tan(e + f*x)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(m + 1)) + 2*b*(c*tan(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**m*hyper((1, m/2 + sympy.S(3)/4), (m/2 + sympy.S(7)/4,), -tan(e + f*x)**2)/(c*f*(2*m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_408():
    f = (d*tan(e + f*x))**m/(a + b*sqrt(c*tan(e + f*x)))
    F = a*(d*tan(e + f*x))**m*(a**2 - b**2*sqrt(-c**2))*tan(e + f*x)*hyper((1, m + 1), (m + 2,), -c*tan(e + f*x)/sqrt(-c**2))/(f*(2*a**4 + 2*b**4*c**2)*(m + 1)) + a*(d*tan(e + f*x))**m*(a**2 + b**2*sqrt(-c**2))*tan(e + f*x)*hyper((1, m + 1), (m + 2,), c*tan(e + f*x)/sqrt(-c**2))/(f*(2*a**4 + 2*b**4*c**2)*(m + 1)) - b*(c*tan(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**m*(a**2 - b**2*sqrt(-c**2))*hyper((1, m + sympy.S(3)/2), (m + sympy.S(5)/2,), -c*tan(e + f*x)/sqrt(-c**2))/(c*f*(a**4 + b**4*c**2)*(2*m + 3)) - b*(c*tan(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**m*(a**2 + b**2*sqrt(-c**2))*hyper((1, m + sympy.S(3)/2), (m + sympy.S(5)/2,), c*tan(e + f*x)/sqrt(-c**2))/(c*f*(a**4 + b**4*c**2)*(2*m + 3)) + b**4*c**2*(d*tan(e + f*x))**m*tan(e + f*x)*hyper((1, 2*m + 2), (2*m + 3,), -b*sqrt(c*tan(e + f*x))/a)/(a*f*(a**4 + b**4*c**2)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_409():
    f = (d*tan(e + f*x))**m/(a + b*sqrt(c*tan(e + f*x)))**2
    F = 4*a**2*b**4*c**2*(d*tan(e + f*x))**m*tan(e + f*x)*hyper((1, 2*m + 2), (2*m + 3,), -b*sqrt(c*tan(e + f*x))/a)/(f*(a**4 + b**4*c**2)**2*(m + 1)) - 2*a*b*(c*tan(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**m*(a**4 - 2*a**2*b**2*sqrt(-c**2) - b**4*c**2)*hyper((1, m + sympy.S(3)/2), (m + sympy.S(5)/2,), -c*tan(e + f*x)/sqrt(-c**2))/(c*f*(a**4 + b**4*c**2)**2*(2*m + 3)) - 2*a*b*(c*tan(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**m*(a**4 + 2*a**2*b**2*sqrt(-c**2) - b**4*c**2)*hyper((1, m + sympy.S(3)/2), (m + sympy.S(5)/2,), c*tan(e + f*x)/sqrt(-c**2))/(c*f*(a**4 + b**4*c**2)**2*(2*m + 3)) + (d*tan(e + f*x))**m*(a**6 - 3*a**4*b**2*sqrt(-c**2) - 3*a**2*b**4*c**2 - b**6*(-c**2)**(sympy.S(3)/2))*tan(e + f*x)*hyper((1, m + 1), (m + 2,), -c*tan(e + f*x)/sqrt(-c**2))/(2*f*(a**4 + b**4*c**2)**2*(m + 1)) + (d*tan(e + f*x))**m*(a**6 + 3*a**4*b**2*sqrt(-c**2) - 3*a**2*b**4*c**2 + b**6*(-c**2)**(sympy.S(3)/2))*tan(e + f*x)*hyper((1, m + 1), (m + 2,), c*tan(e + f*x)/sqrt(-c**2))/(2*f*(a**4 + b**4*c**2)**2*(m + 1)) + b**4*c**2*(d*tan(e + f*x))**m*tan(e + f*x)*hyper((2, 2*m + 2), (2*m + 3,), -b*sqrt(c*tan(e + f*x))/a)/(a**2*f*(a**4 + b**4*c**2)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_410():
    f = (b*(c*tan(e + f*x))**n)**p*(d*tan(e + f*x))**m
    F = (b*(c*tan(e + f*x))**n)**p*(d*tan(e + f*x))**m*tan(e + f*x)*hyper((1, m/2 + n*p/2 + sympy.S.Half), (m/2 + n*p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(m + n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_411():
    f = (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)**2
    F = (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)**3*hyper((1, n*p/2 + sympy.S(3)/2), (n*p/2 + sympy.S(5)/2,), -tan(e + f*x)**2)/(f*(n*p + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_412():
    f = (b*(c*tan(e + f*x))**n)**p
    F = (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)*hyper((1, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_413():
    f = (b*(c*tan(e + f*x))**n)**p*cot(e + f*x)**2
    F = -(b*(c*tan(e + f*x))**n)**p*cot(e + f*x)*hyper((1, n*p/2 + sympy.S(-1)/2), (n*p/2 + sympy.S.Half,), -tan(e + f*x)**2)/(f*(-n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_414():
    f = (b*(c*tan(e + f*x))**n)**p*cot(e + f*x)**4
    F = -(b*(c*tan(e + f*x))**n)**p*cot(e + f*x)**3*hyper((1, n*p/2 + sympy.S(-3)/2), (n*p/2 + sympy.S(-1)/2,), -tan(e + f*x)**2)/(f*(-n*p + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_415():
    f = (b*(c*tan(e + f*x))**n)**p*cot(e + f*x)**6
    F = -(b*(c*tan(e + f*x))**n)**p*cot(e + f*x)**5*hyper((1, n*p/2 + sympy.S(-5)/2), (n*p/2 + sympy.S(-3)/2,), -tan(e + f*x)**2)/(f*(-n*p + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_416():
    f = (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)**3
    F = (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)**4*hyper((1, n*p/2 + 2), (n*p/2 + 3,), -tan(e + f*x)**2)/(f*(n*p + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_417():
    f = (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)
    F = (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)**2*hyper((1, n*p/2 + 1), (n*p/2 + 2,), -tan(e + f*x)**2)/(f*(n*p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_418():
    f = (b*(c*tan(e + f*x))**n)**p*cot(e + f*x)
    F = (b*(c*tan(e + f*x))**n)**p*hyper((1, n*p/2), (n*p/2 + 1,), -tan(e + f*x)**2)/(f*n*p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_419():
    f = (b*(c*tan(e + f*x))**n)**p*cot(e + f*x)**3
    F = -(b*(c*tan(e + f*x))**n)**p*cot(e + f*x)**2*hyper((1, n*p/2 - 1), (n*p/2,), -tan(e + f*x)**2)/(f*(-n*p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_420():
    f = (d*tan(e + f*x))**m*(a + b*(c*tan(e + f*x))**n)**p
    F = sympy.Function('Unintegrable')((((Symbol('d') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_421():
    f = (b*tan(e + f*x)**2)**p*(d*cot(e + f*x))**m
    F = (b*tan(e + f*x)**2)**p*(d*cot(e + f*x))**m*tan(e + f*x)*hyper((1, -m/2 + p + sympy.S.Half), (-m/2 + p + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(-m + 2*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_422():
    f = (d*cot(e + f*x))**m*(a + b*tan(e + f*x)**2)**p
    F = (d*cot(e + f*x))**m*(a + b*tan(e + f*x)**2)**p*tan(e + f*x)*appellf1(sympy.S.Half - m/2, 1, -p, sympy.S(3)/2 - m/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/a)/(f*(1 - m)*(1 + b*tan(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_423():
    f = (b*(c*tan(e + f*x))**n)**p*(d*cot(e + f*x))**m
    F = (b*(c*tan(e + f*x))**n)**p*(d*cot(e + f*x))**m*tan(e + f*x)*hyper((1, -m/2 + n*p/2 + sympy.S.Half), (-m/2 + n*p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(-m + n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_424():
    f = (d*cot(e + f*x))**m*(a + b*(c*tan(e + f*x))**n)**p
    F = ((Symbol('d') * sympy.cot((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * ((sympy.tan((Symbol('e') + (Symbol('f') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')) * sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p')) * (((sympy.tan((Symbol('e') + (Symbol('f') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_425():
    f = (a + b*tan(c + d*x)**2)*sec(c + d*x)**3
    F = b*tan(c + d*x)*sec(c + d*x)**3/(4*d) + (4*a - b)*tan(c + d*x)*sec(c + d*x)/(8*d) + (4*a - b)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_426():
    f = (a + b*tan(c + d*x)**2)*sec(c + d*x)
    F = b*tan(c + d*x)*sec(c + d*x)/(2*d) + (2*a - b)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_427():
    f = (a + b*tan(c + d*x)**2)*cos(c + d*x)
    F = b*atanh(sin(c + d*x))/d + (a - b)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_428():
    f = (a + b*tan(c + d*x)**2)*cos(c + d*x)**3
    F = a*sin(c + d*x)/d - (a - b)*sin(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_429():
    f = (a + b*tan(c + d*x)**2)*cos(c + d*x)**5
    F = a*sin(c + d*x)/d + (a - b)*sin(c + d*x)**5/(5*d) - (2*a - b)*sin(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_430():
    f = (a + b*tan(c + d*x)**2)*cos(c + d*x)**7
    F = a*sin(c + d*x)/d - (a - b)*sin(c + d*x)**7/(7*d) + (3*a - 2*b)*sin(c + d*x)**5/(5*d) - (3*a - b)*sin(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_431():
    f = (a + b*tan(c + d*x)**2)*sec(c + d*x)**6
    F = a*tan(c + d*x)/d + b*tan(c + d*x)**7/(7*d) + (a + 2*b)*tan(c + d*x)**5/(5*d) + (2*a + b)*tan(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_432():
    f = (a + b*tan(c + d*x)**2)*sec(c + d*x)**4
    F = a*tan(c + d*x)/d + b*tan(c + d*x)**5/(5*d) + (a + b)*tan(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_433():
    f = (a + b*tan(c + d*x)**2)*sec(c + d*x)**2
    F = a*tan(c + d*x)/d + b*tan(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_434():
    f = (a + b*tan(c + d*x)**2)*cos(c + d*x)**2
    F = x*(a/2 + b/2) + (a - b)*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_435():
    f = (a + b*tan(c + d*x)**2)*cos(c + d*x)**4
    F = x*(3*a/8 + b/8) + (a - b)*sin(c + d*x)*cos(c + d*x)**3/(4*d) + (3*a + b)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_436():
    f = (a + b*tan(c + d*x)**2)*cos(c + d*x)**6
    F = x*(5*a/16 + b/16) + (a - b)*sin(c + d*x)*cos(c + d*x)**5/(6*d) + (5*a + b)*sin(c + d*x)*cos(c + d*x)**3/(24*d) + (5*a + b)*sin(c + d*x)*cos(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_437():
    f = (a + b*tan(c + d*x)**2)**2*sec(c + d*x)**3
    F = b*(a - (a - b)*sin(c + d*x)**2)*tan(c + d*x)*sec(c + d*x)**5/(6*d) + b*(8*a - 3*b)*tan(c + d*x)*sec(c + d*x)**3/(24*d) + (8*a**2 - 4*a*b + b**2)*tan(c + d*x)*sec(c + d*x)/(16*d) + (8*a**2 - 4*a*b + b**2)*atanh(sin(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_438():
    f = (a + b*tan(c + d*x)**2)**2*sec(c + d*x)
    F = b*(a - (a - b)*sin(c + d*x)**2)*tan(c + d*x)*sec(c + d*x)**3/(4*d) + b*(6*a - 3*b)*tan(c + d*x)*sec(c + d*x)/(8*d) + (8*a**2 - 8*a*b + 3*b**2)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_439():
    f = (a + b*tan(c + d*x)**2)**2*cos(c + d*x)
    F = b**2*tan(c + d*x)*sec(c + d*x)/(2*d) + b*(4*a - 3*b)*atanh(sin(c + d*x))/(2*d) + (a - b)**2*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_440():
    f = (a + b*tan(c + d*x)**2)**2*cos(c + d*x)**3
    F = b**2*atanh(sin(c + d*x))/d - (a - b)**2*sin(c + d*x)**3/(3*d) + (a**2 - b**2)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_441():
    f = (a + b*tan(c + d*x)**2)**2*cos(c + d*x)**5
    F = a**2*sin(c + d*x)/d - 2*a*(a - b)*sin(c + d*x)**3/(3*d) + (a - b)**2*sin(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_442():
    f = (a + b*tan(c + d*x)**2)**2*cos(c + d*x)**7
    F = a**2*sin(c + d*x)/d - a*(3*a - 2*b)*sin(c + d*x)**3/(3*d) - (a - b)**2*sin(c + d*x)**7/(7*d) + (a - b)*(3*a - b)*sin(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_443():
    f = (a + b*tan(c + d*x)**2)**2*cos(c + d*x)**9
    F = a**2*sin(c + d*x)/d - 2*a*(2*a - b)*sin(c + d*x)**3/(3*d) + (a - b)**2*sin(c + d*x)**9/(9*d) - (2*a - 2*b)*(2*a - b)*sin(c + d*x)**7/(7*d) + (6*a**2 - 6*a*b + b**2)*sin(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_444():
    f = (a + b*tan(c + d*x)**2)**2*sec(c + d*x)**6
    F = a**2*tan(c + d*x)/d + 2*a*(a + b)*tan(c + d*x)**3/(3*d) + b**2*tan(c + d*x)**9/(9*d) + 2*b*(a + b)*tan(c + d*x)**7/(7*d) + (a**2 + 4*a*b + b**2)*tan(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_445():
    f = (a + b*tan(c + d*x)**2)**2*sec(c + d*x)**4
    F = a**2*tan(c + d*x)/d + a*(a + 2*b)*tan(c + d*x)**3/(3*d) + b**2*tan(c + d*x)**7/(7*d) + b*(2*a + b)*tan(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_446():
    f = (a + b*tan(c + d*x)**2)**2*sec(c + d*x)**2
    F = a**2*tan(c + d*x)/d + 2*a*b*tan(c + d*x)**3/(3*d) + b**2*tan(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_447():
    f = (a + b*tan(c + d*x)**2)**2*cos(c + d*x)**2
    F = b**2*tan(c + d*x)/d + x*(a/2 - b/2)*(a + 3*b) + (a - b)**2*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_448():
    f = (a + b*tan(c + d*x)**2)**2*cos(c + d*x)**4
    F = x*(3*a**2/8 + a*b/4 + 3*b**2/8) + (a - b)*(a + b*tan(c + d*x)**2)*sin(c + d*x)*cos(c + d*x)**3/(4*d) + (3*a**2 - 3*b**2)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_449():
    f = (a + b*tan(c + d*x)**2)**2*cos(c + d*x)**6
    F = x*(5*a**2/16 + a*b/8 + b**2/16) + (a - b)*(a + b*tan(c + d*x)**2)*sin(c + d*x)*cos(c + d*x)**5/(6*d) + (a - b)*(5*a + 3*b)*sin(c + d*x)*cos(c + d*x)**3/(24*d) + (5*a**2 + 2*a*b + b**2)*sin(c + d*x)*cos(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_450():
    f = sec(c + d*x)**5/(a + b*tan(c + d*x)**2)
    F = tan(c + d*x)*sec(c + d*x)/(2*b*d) - (2*a - 3*b)*atanh(sin(c + d*x))/(2*b**2*d) + (a - b)**(sympy.S(3)/2)*atanh(sqrt(a - b)*sin(c + d*x)/sqrt(a))/(sqrt(a)*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_451():
    f = sec(c + d*x)**3/(a + b*tan(c + d*x)**2)
    F = atanh(sin(c + d*x))/(b*d) - sqrt(a - b)*atanh(sqrt(a - b)*sin(c + d*x)/sqrt(a))/(sqrt(a)*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_452():
    f = sec(c + d*x)/(a + b*tan(c + d*x)**2)
    F = atanh(sqrt(a - b)*sin(c + d*x)/sqrt(a))/(sqrt(a)*d*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_453():
    f = cos(c + d*x)/(a + b*tan(c + d*x)**2)
    F = sin(c + d*x)/(d*(a - b)) - b*atanh(sqrt(a - b)*sin(c + d*x)/sqrt(a))/(sqrt(a)*d*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_454():
    f = cos(c + d*x)**3/(a + b*tan(c + d*x)**2)
    F = (a - 2*b)*sin(c + d*x)/(d*(a - b)**2) - sin(c + d*x)**3/(d*(3*a - 3*b)) + b**2*atanh(sqrt(a - b)*sin(c + d*x)/sqrt(a))/(sqrt(a)*d*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_455():
    f = cos(c + d*x)**5/(a + b*tan(c + d*x)**2)
    F = sin(c + d*x)**5/(d*(5*a - 5*b)) - (2*a - 3*b)*sin(c + d*x)**3/(3*d*(a - b)**2) + (a**2 - 3*a*b + 3*b**2)*sin(c + d*x)/(d*(a - b)**3) - b**3*atanh(sqrt(a - b)*sin(c + d*x)/sqrt(a))/(sqrt(a)*d*(a - b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_456():
    f = sec(c + d*x)**8/(a + b*tan(c + d*x)**2)
    F = tan(c + d*x)**5/(5*b*d) - (a - 3*b)*tan(c + d*x)**3/(3*b**2*d) + (a**2 - 3*a*b + 3*b**2)*tan(c + d*x)/(b**3*d) - (a - b)**3*atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(sqrt(a)*b**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_457():
    f = sec(c + d*x)**6/(a + b*tan(c + d*x)**2)
    F = tan(c + d*x)**3/(3*b*d) - (a - 2*b)*tan(c + d*x)/(b**2*d) + (a - b)**2*atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(sqrt(a)*b**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_458():
    f = sec(c + d*x)**4/(a + b*tan(c + d*x)**2)
    F = tan(c + d*x)/(b*d) - (a - b)*atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(sqrt(a)*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_459():
    f = sec(c + d*x)**2/(a + b*tan(c + d*x)**2)
    F = atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(sqrt(a)*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_460():
    f = cos(c + d*x)**2/(a + b*tan(c + d*x)**2)
    F = x*(a - 3*b)/(2*(a - b)**2) + sin(c + d*x)*cos(c + d*x)/(d*(2*a - 2*b)) + b**(sympy.S(3)/2)*atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(sqrt(a)*d*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_461():
    f = cos(c + d*x)**4/(a + b*tan(c + d*x)**2)
    F = x*(3*a**2 - 10*a*b + 15*b**2)/(8*(a - b)**3) + sin(c + d*x)*cos(c + d*x)**3/(d*(4*a - 4*b)) + (3*a - 7*b)*sin(c + d*x)*cos(c + d*x)/(8*d*(a - b)**2) - b**(sympy.S(5)/2)*atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(sqrt(a)*d*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_462():
    f = sec(c + d*x)**7/(a + b*tan(c + d*x)**2)**2
    F = tan(c + d*x)*sec(c + d*x)/(2*b*d*(a - (a - b)*sin(c + d*x)**2)) - (4*a - 5*b)*atanh(sin(c + d*x))/(2*b**3*d) + (a - b)*(2*a - b)*sin(c + d*x)/(2*a*b**2*d*(a - (a - b)*sin(c + d*x)**2)) + (a - b)**(sympy.S(3)/2)*(4*a + b)*atanh(sqrt(a - b)*sin(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_463():
    f = sec(c + d*x)**5/(a + b*tan(c + d*x)**2)**2
    F = atanh(sin(c + d*x))/(b**2*d) - (a - b)*sin(c + d*x)/(2*a*b*d*(a - (a - b)*sin(c + d*x)**2)) - sqrt(a - b)*(2*a + b)*atanh(sqrt(a - b)*sin(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_464():
    f = sec(c + d*x)**3/(a + b*tan(c + d*x)**2)**2
    F = sin(c + d*x)/(2*a*d*(a - (a - b)*sin(c + d*x)**2)) + atanh(sqrt(a - b)*sin(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_465():
    f = sec(c + d*x)/(a + b*tan(c + d*x)**2)**2
    F = -b*sin(c + d*x)/(2*a*d*(a - b)*(a - (a - b)*sin(c + d*x)**2)) + (2*a - b)*atanh(sqrt(a - b)*sin(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_466():
    f = cos(c + d*x)/(a + b*tan(c + d*x)**2)**2
    F = sin(c + d*x)/(d*(a - b)**2) + b**2*sin(c + d*x)/(2*a*d*(a - b)**2*(a - (a - b)*sin(c + d*x)**2)) - b*(4*a - b)*atanh(sqrt(a - b)*sin(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_467():
    f = cos(c + d*x)**3/(a + b*tan(c + d*x)**2)**2
    F = (a - 3*b)*sin(c + d*x)/(d*(a - b)**3) - sin(c + d*x)**3/(3*d*(a - b)**2) - b**3*sin(c + d*x)/(2*a*d*(a - b)**3*(a - (a - b)*sin(c + d*x)**2)) + b**2*(6*a - b)*atanh(sqrt(a - b)*sin(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a - b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_468():
    f = sec(c + d*x)**8/(a + b*tan(c + d*x)**2)**2
    F = tan(c + d*x)**3/(3*b**2*d) - (2*a - 3*b)*tan(c + d*x)/(b**3*d) - (a - b)**3*tan(c + d*x)/(2*a*b**3*d*(a + b*tan(c + d*x)**2)) + (a - b)**2*(5*a + b)*atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_469():
    f = sec(c + d*x)**6/(a + b*tan(c + d*x)**2)**2
    F = tan(c + d*x)/(b**2*d) + (a - b)**2*tan(c + d*x)/(2*a*b**2*d*(a + b*tan(c + d*x)**2)) - (3*a**2 - 2*a*b - b**2)*atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_470():
    f = sec(c + d*x)**4/(a + b*tan(c + d*x)**2)**2
    F = -(a - b)*tan(c + d*x)/(2*a*b*d*(a + b*tan(c + d*x)**2)) + (a + b)*atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_471():
    f = sec(c + d*x)**2/(a + b*tan(c + d*x)**2)**2
    F = tan(c + d*x)/(2*a*d*(a + b*tan(c + d*x)**2)) + atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_472():
    f = cos(c + d*x)**2/(a + b*tan(c + d*x)**2)**2
    F = x*(a - 5*b)/(2*(a - b)**3) + sin(c + d*x)*cos(c + d*x)/(d*(a + b*tan(c + d*x)**2)*(2*a - 2*b)) + b*(a + b)*tan(c + d*x)/(2*a*d*(a - b)**2*(a + b*tan(c + d*x)**2)) + b**(sympy.S(3)/2)*(5*a - b)*atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_473():
    f = cos(c + d*x)**4/(a + b*tan(c + d*x)**2)**2
    F = x*(3*a**2 - 14*a*b + 35*b**2)/(8*(a - b)**4) + sin(c + d*x)*cos(c + d*x)**3/(d*(a + b*tan(c + d*x)**2)*(4*a - 4*b)) + (3*a - 9*b)*sin(c + d*x)*cos(c + d*x)/(8*d*(a - b)**2*(a + b*tan(c + d*x)**2)) + b*(a - 4*b)*(3*a + b)*tan(c + d*x)/(8*a*d*(a - b)**3*(a + b*tan(c + d*x)**2)) - b**(sympy.S(5)/2)*(7*a - b)*atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a - b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_474():
    f = (b*tan(e + f*x)**2)**p*(d*sec(e + f*x))**m
    F = (b*tan(e + f*x)**2)**p*(d*sec(e + f*x))**m*(cos(e + f*x)**2)**(m/2 + p + sympy.S.Half)*tan(e + f*x)*hyper((p + sympy.S.Half, m/2 + p + sympy.S.Half), (p + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_475():
    f = (d*sec(e + f*x))**m*(a + b*tan(e + f*x)**2)**p
    F = (d*sec(e + f*x))**m*(a + b*tan(e + f*x)**2)**p*tan(e + f*x)*appellf1(sympy.S.Half, -p, 1 - m/2, sympy.S(3)/2, -b*tan(e + f*x)**2/a, -tan(e + f*x)**2)/(f*(1 + b*tan(e + f*x)**2/a)**p*(sec(e + f*x)**2)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_476():
    f = (b*(c*tan(e + f*x))**n)**p*(d*sec(e + f*x))**m
    F = (b*(c*tan(e + f*x))**n)**p*(d*sec(e + f*x))**m*(cos(e + f*x)**2)**(m/2 + n*p/2 + sympy.S.Half)*tan(e + f*x)*hyper((n*p/2 + sympy.S.Half, m/2 + n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_477():
    f = (b*(c*tan(e + f*x))**n)**p*sec(e + f*x)**6
    F = (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)**5/(f*(n*p + 5)) + 2*(b*(c*tan(e + f*x))**n)**p*tan(e + f*x)**3/(f*(n*p + 3)) + (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_478():
    f = (b*(c*tan(e + f*x))**n)**p*sec(e + f*x)**4
    F = (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)**3/(f*(n*p + 3)) + (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_479():
    f = (b*(c*tan(e + f*x))**n)**p*sec(e + f*x)**2
    F = (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_480():
    f = (b*(c*tan(e + f*x))**n)**p
    F = (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)*hyper((1, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_481():
    f = (b*(c*tan(e + f*x))**n)**p*cos(e + f*x)**2
    F = (b*(c*tan(e + f*x))**n)**p*tan(e + f*x)*hyper((2, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_482():
    f = (b*(c*tan(e + f*x))**n)**p*sec(e + f*x)**3
    F = (b*(c*tan(e + f*x))**n)**p*(cos(e + f*x)**2)**(n*p/2 + 2)*tan(e + f*x)*hyper((n*p/2 + sympy.S.Half, n*p/2 + 2), (n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)*sec(e + f*x)**3/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_483():
    f = (b*(c*tan(e + f*x))**n)**p*sec(e + f*x)
    F = (b*(c*tan(e + f*x))**n)**p*(cos(e + f*x)**2)**(n*p/2 + 1)*tan(e + f*x)*hyper((n*p/2 + sympy.S.Half, n*p/2 + 1), (n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)*sec(e + f*x)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_484():
    f = (b*(c*tan(e + f*x))**n)**p*cos(e + f*x)
    F = (b*(c*tan(e + f*x))**n)**p*(cos(e + f*x)**2)**(n*p/2)*sin(e + f*x)*hyper((n*p/2, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_485():
    f = (b*(c*tan(e + f*x))**n)**p*cos(e + f*x)**3
    F = (b*(c*tan(e + f*x))**n)**p*(cos(e + f*x)**2)**(n*p/2)*sin(e + f*x)*hyper((n*p/2 - 1, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_486():
    f = (d*sec(e + f*x))**m*(a + b*(c*tan(e + f*x))**n)**p
    F = sympy.Function('Unintegrable')((((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_487():
    f = (a + b*(c*tan(e + f*x))**n)**p*sec(e + f*x)**3
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_488():
    f = (a + b*(c*tan(e + f*x))**n)**p*sec(e + f*x)
    F = sympy.Function('Unintegrable')((sympy.sec((Symbol('e') + (Symbol('f') * x))) * ((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_489():
    f = (a + b*(c*tan(e + f*x))**n)**p*cos(e + f*x)
    F = sympy.Function('Unintegrable')((sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_490():
    f = (a + b*(c*tan(e + f*x))**n)**p*cos(e + f*x)**3
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_491():
    f = (a + b*(c*tan(e + f*x))**n)**p*sec(e + f*x)**6
    F = (a + b*(c*tan(e + f*x))**n)**p*tan(e + f*x)**5*hyper((-p, 5/n), ((n + 5)/n,), -b*(c*tan(e + f*x))**n/a)/(5*f*(1 + b*(c*tan(e + f*x))**n/a)**p) + 2*(a + b*(c*tan(e + f*x))**n)**p*tan(e + f*x)**3*hyper((-p, 3/n), ((n + 3)/n,), -b*(c*tan(e + f*x))**n/a)/(3*f*(1 + b*(c*tan(e + f*x))**n/a)**p) + (a + b*(c*tan(e + f*x))**n)**p*tan(e + f*x)*hyper((1/n, -p), (1 + 1/n,), -b*(c*tan(e + f*x))**n/a)/(f*(1 + b*(c*tan(e + f*x))**n/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_492():
    f = (a + b*(c*tan(e + f*x))**n)**p*sec(e + f*x)**4
    F = (a + b*(c*tan(e + f*x))**n)**p*tan(e + f*x)**3*hyper((-p, 3/n), ((n + 3)/n,), -b*(c*tan(e + f*x))**n/a)/(3*f*(1 + b*(c*tan(e + f*x))**n/a)**p) + (a + b*(c*tan(e + f*x))**n)**p*tan(e + f*x)*hyper((1/n, -p), (1 + 1/n,), -b*(c*tan(e + f*x))**n/a)/(f*(1 + b*(c*tan(e + f*x))**n/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_493():
    f = (a + b*(c*tan(e + f*x))**n)**p*sec(e + f*x)**2
    F = (a + b*(c*tan(e + f*x))**n)**p*tan(e + f*x)*hyper((1/n, -p), (1 + 1/n,), -b*(c*tan(e + f*x))**n/a)/(f*(1 + b*(c*tan(e + f*x))**n/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_494():
    f = (a + b*(c*tan(e + f*x))**n)**p
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p')), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_495():
    f = (a + b*(c*tan(e + f*x))**n)**p*cos(e + f*x)**2
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_496():
    f = (b*tan(e + f*x)**2)**p*(d*csc(e + f*x))**m
    F = (b*tan(e + f*x)**2)**p*(d*csc(e + f*x))**m*(cos(e + f*x)**2)**(p + sympy.S.Half)*tan(e + f*x)*hyper((p + sympy.S.Half, -m/2 + p + sympy.S.Half), (-m/2 + p + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(-m + 2*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_497():
    f = (d*csc(e + f*x))**m*(a + b*tan(e + f*x)**2)**p
    F = (d*csc(e + f*x))**m*(a + b*tan(e + f*x)**2)**p*tan(e + f*x)*appellf1(sympy.S.Half - m/2, -p, 1 - m/2, sympy.S(3)/2 - m/2, -b*tan(e + f*x)**2/a, -tan(e + f*x)**2)/(f*(1 - m)*(1 + b*tan(e + f*x)**2/a)**p*(sec(e + f*x)**2)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_498():
    f = (b*(c*tan(e + f*x))**n)**p*(d*csc(e + f*x))**m
    F = (b*(c*tan(e + f*x))**n)**p*(d*csc(e + f*x))**m*(cos(e + f*x)**2)**(n*p/2 + sympy.S.Half)*tan(e + f*x)*hyper((n*p/2 + sympy.S.Half, -m/2 + n*p/2 + sympy.S.Half), (-m/2 + n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(-m + n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_7_d_trig_pow_m_a_plus_b_c_tan_pow_n_pow_p_499():
    f = (d*csc(e + f*x))**m*(a + b*(c*tan(e + f*x))**n)**p
    F = ((Symbol('d') * sympy.csc((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * ((sympy.sin((Symbol('e') + (Symbol('f') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')) * sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p')) * (((sympy.sin((Symbol('e') + (Symbol('f') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')))**(Integer(-1))), x)
    assert integrate(f, x) == F

