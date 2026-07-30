"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.4 Cotangent/4.4.7 (d trig)^m (a+b (c cot)^n)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

A, C, a, b, c, d = symbols('A C a b c d')

def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_1():
    f = (A + C*cot(c + d*x)**2)/sqrt(b*tan(c + d*x))
    F = -2*C*b/(3*d*(b*tan(c + d*x))**(sympy.S(3)/2)) - sqrt(2)*(A - C)*log(sqrt(b)*tan(c + d*x) + sqrt(b) - sqrt(2)*sqrt(b*tan(c + d*x)))/(4*sqrt(b)*d) + sqrt(2)*(A - C)*log(sqrt(b)*tan(c + d*x) + sqrt(b) + sqrt(2)*sqrt(b*tan(c + d*x)))/(4*sqrt(b)*d) - sqrt(2)*(A - C)*atan(1 - sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*sqrt(b)*d) + sqrt(2)*(A - C)*atan(1 + sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_2():
    f = a + b*cot(c + d*x)**2
    F = a*x - b*x - b*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_3():
    f = (a + b*cot(c + d*x)**2)**2
    F = -b**2*cot(c + d*x)**3/(3*d) - b*(2*a - b)*cot(c + d*x)/d + x*(a - b)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_4():
    f = (a + b*cot(c + d*x)**2)**3
    F = -b**3*cot(c + d*x)**5/(5*d) - b**2*(3*a - b)*cot(c + d*x)**3/(3*d) - b*(3*a**2 - 3*a*b + b**2)*cot(c + d*x)/d + x*(a - b)**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_5():
    f = 1/(a + b*cot(c + d*x)**2)
    F = x/(a - b) + sqrt(b)*atan(sqrt(b)*cot(c + d*x)/sqrt(a))/(sqrt(a)*d*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_6():
    f = (a + b*cot(c + d*x)**2)**(-2)
    F = x/(a - b)**2 + b*cot(c + d*x)/(2*a*d*(a - b)*(a + b*cot(c + d*x)**2)) + sqrt(b)*(3*a - b)*atan(sqrt(b)*cot(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_7():
    f = (a + b*cot(c + d*x)**2)**(-3)
    F = x/(a - b)**3 + b*cot(c + d*x)/(4*a*d*(a - b)*(a + b*cot(c + d*x)**2)**2) + b*(7*a - 3*b)*cot(c + d*x)/(8*a**2*d*(a - b)**2*(a + b*cot(c + d*x)**2)) + sqrt(b)*(15*a**2 - 10*a*b + 3*b**2)*atan(sqrt(b)*cot(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_8():
    f = (cot(x)**2 + 1)**(sympy.S(3)/2)
    F = -sqrt(csc(x)**2)*cot(x)/2 - asinh(cot(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_9():
    f = sqrt(cot(x)**2 + 1)
    F = -asinh(cot(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_10():
    f = 1/sqrt(cot(x)**2 + 1)
    F = -cot(x)/sqrt(csc(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_11():
    f = (-cot(x)**2 - 1)**(sympy.S(3)/2)
    F = sqrt(-csc(x)**2)*cot(x)/2 - atan(cot(x)/sqrt(-csc(x)**2))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_12():
    f = sqrt(-cot(x)**2 - 1)
    F = atan(cot(x)/sqrt(-csc(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_13():
    f = 1/sqrt(-cot(x)**2 - 1)
    F = -cot(x)/sqrt(-csc(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_14():
    f = cot(x)**3/sqrt(a*cot(x)**2 + a)
    F = -1/sqrt(a*csc(x)**2) - sqrt(a*csc(x)**2)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_15():
    f = cot(x)**2/sqrt(a*cot(x)**2 + a)
    F = cot(x)/sqrt(a*csc(x)**2) - atanh(cos(x))*csc(x)/sqrt(a*csc(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_16():
    f = cot(x)/sqrt(a*cot(x)**2 + a)
    F = 1/sqrt(a*csc(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_17():
    f = tan(x)/sqrt(a*cot(x)**2 + a)
    F = -1/sqrt(a*csc(x)**2) + atanh(sqrt(a*csc(x)**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_18():
    f = tan(x)**2/sqrt(a*cot(x)**2 + a)
    F = cot(x)/sqrt(a*csc(x)**2) + csc(x)*sec(x)/sqrt(a*csc(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_19():
    f = sqrt(a + b*cot(x)**2)*cot(x)**3
    F = -sqrt(a - b)*atanh(sqrt(a + b*cot(x)**2)/sqrt(a - b)) + sqrt(a + b*cot(x)**2) - (a + b*cot(x)**2)**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_20():
    f = sqrt(a + b*cot(x)**2)*cot(x)
    F = sqrt(a - b)*atanh(sqrt(a + b*cot(x)**2)/sqrt(a - b)) - sqrt(a + b*cot(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_21():
    f = sqrt(a + b*cot(x)**2)*tan(x)
    F = sqrt(a)*atanh(sqrt(a + b*cot(x)**2)/sqrt(a)) - sqrt(a - b)*atanh(sqrt(a + b*cot(x)**2)/sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_22():
    f = sqrt(a + b*cot(x)**2)*cot(x)**2
    F = sqrt(a - b)*atan(sqrt(a - b)*cot(x)/sqrt(a + b*cot(x)**2)) - sqrt(a + b*cot(x)**2)*cot(x)/2 - (a - 2*b)*atanh(sqrt(b)*cot(x)/sqrt(a + b*cot(x)**2))/(2*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_23():
    f = sqrt(a + b*cot(x)**2)
    F = -sqrt(b)*atanh(sqrt(b)*cot(x)/sqrt(a + b*cot(x)**2)) - sqrt(a - b)*atan(sqrt(a - b)*cot(x)/sqrt(a + b*cot(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_24():
    f = sqrt(a + b*cot(x)**2)*tan(x)**2
    F = sqrt(a - b)*atan(sqrt(a - b)*cot(x)/sqrt(a + b*cot(x)**2)) + sqrt(a + b*cot(x)**2)*tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_25():
    f = sqrt(a + b*cot(x)**2)*tan(x)**4
    F = -sqrt(a - b)*atan(sqrt(a - b)*cot(x)/sqrt(a + b*cot(x)**2)) + sqrt(a + b*cot(x)**2)*tan(x)**3/3 - sqrt(a + b*cot(x)**2)*(3*a - b)*tan(x)/(3*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_26():
    f = (a + b*cot(x)**2)**(sympy.S(3)/2)*cot(x)**3
    F = -(a - b)**(sympy.S(3)/2)*atanh(sqrt(a + b*cot(x)**2)/sqrt(a - b)) + (a - b)*sqrt(a + b*cot(x)**2) + (a + b*cot(x)**2)**(sympy.S(3)/2)/3 - (a + b*cot(x)**2)**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_27():
    f = (a + b*cot(x)**2)**(sympy.S(3)/2)*cot(x)**2
    F = -b*sqrt(a + b*cot(x)**2)*cot(x)**3/4 - (5*a/8 - b/2)*sqrt(a + b*cot(x)**2)*cot(x) + (a - b)**(sympy.S(3)/2)*atan(sqrt(a - b)*cot(x)/sqrt(a + b*cot(x)**2)) - (3*a**2 - 12*a*b + 8*b**2)*atanh(sqrt(b)*cot(x)/sqrt(a + b*cot(x)**2))/(8*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_28():
    f = (a + b*cot(x)**2)**(sympy.S(3)/2)*cot(x)
    F = (a - b)**(sympy.S(3)/2)*atanh(sqrt(a + b*cot(x)**2)/sqrt(a - b)) - (a - b)*sqrt(a + b*cot(x)**2) - (a + b*cot(x)**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_29():
    f = (a + b*cot(x)**2)**(sympy.S(3)/2)*tan(x)
    F = a**(sympy.S(3)/2)*atanh(sqrt(a + b*cot(x)**2)/sqrt(a)) - b*sqrt(a + b*cot(x)**2) - (a - b)**(sympy.S(3)/2)*atanh(sqrt(a + b*cot(x)**2)/sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_30():
    f = (a + b*cot(x)**2)**(sympy.S(3)/2)*tan(x)**2
    F = a*sqrt(a + b*cot(x)**2)*tan(x) - b**(sympy.S(3)/2)*atanh(sqrt(b)*cot(x)/sqrt(a + b*cot(x)**2)) + (a - b)**(sympy.S(3)/2)*atan(sqrt(a - b)*cot(x)/sqrt(a + b*cot(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_31():
    f = (a + b*cot(c + d*x)**2)**(sympy.S(5)/2)
    F = -sqrt(b)*(15*a**2 - 20*a*b + 8*b**2)*atanh(sqrt(b)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2))/(8*d) - b*(a + b*cot(c + d*x)**2)**(sympy.S(3)/2)*cot(c + d*x)/(4*d) - b*sqrt(a + b*cot(c + d*x)**2)*(7*a - 4*b)*cot(c + d*x)/(8*d) - (a - b)**(sympy.S(5)/2)*atan(sqrt(a - b)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_32():
    f = (a + b*cot(c + d*x)**2)**(sympy.S(3)/2)
    F = -sqrt(b)*(3*a - 2*b)*atanh(sqrt(b)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2))/(2*d) - b*sqrt(a + b*cot(c + d*x)**2)*cot(c + d*x)/(2*d) - (a - b)**(sympy.S(3)/2)*atan(sqrt(a - b)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_33():
    f = sqrt(a + b*cot(c + d*x)**2)
    F = -sqrt(b)*atanh(sqrt(b)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2))/d - sqrt(a - b)*atan(sqrt(a - b)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_34():
    f = 1/sqrt(a + b*cot(c + d*x)**2)
    F = -atan(sqrt(a - b)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2))/(d*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_35():
    f = (a + b*cot(c + d*x)**2)**(sympy.S(-3)/2)
    F = -atan(sqrt(a - b)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2))/(d*(a - b)**(sympy.S(3)/2)) + b*cot(c + d*x)/(a*d*(a - b)*sqrt(a + b*cot(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_36():
    f = (a + b*cot(c + d*x)**2)**(sympy.S(-5)/2)
    F = -atan(sqrt(a - b)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2))/(d*(a - b)**(sympy.S(5)/2)) + b*cot(c + d*x)/(3*a*d*(a - b)*(a + b*cot(c + d*x)**2)**(sympy.S(3)/2)) + b*(5*a - 2*b)*cot(c + d*x)/(3*a**2*d*(a - b)**2*sqrt(a + b*cot(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_37():
    f = (a + b*cot(c + d*x)**2)**(sympy.S(-7)/2)
    F = -atan(sqrt(a - b)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2))/(d*(a - b)**(sympy.S(7)/2)) + b*cot(c + d*x)/(5*a*d*(a - b)*(a + b*cot(c + d*x)**2)**(sympy.S(5)/2)) + b*(9*a - 4*b)*cot(c + d*x)/(15*a**2*d*(a - b)**2*(a + b*cot(c + d*x)**2)**(sympy.S(3)/2)) + b*(33*a**2 - 26*a*b + 8*b**2)*cot(c + d*x)/(15*a**3*d*(a - b)**3*sqrt(a + b*cot(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_38():
    f = (1 - cot(x)**2)**(sympy.S(3)/2)
    F = sqrt(1 - cot(x)**2)*cot(x)/2 + 5*asin(cot(x))/2 - 2*sqrt(2)*atan(sqrt(2)*cot(x)/sqrt(1 - cot(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_39():
    f = sqrt(1 - cot(x)**2)
    F = asin(cot(x)) - sqrt(2)*atan(sqrt(2)*cot(x)/sqrt(1 - cot(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_40():
    f = 1/sqrt(1 - cot(x)**2)
    F = -sqrt(2)*atan(sqrt(2)*cot(x)/sqrt(1 - cot(x)**2))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_41():
    f = (cot(x)**2 - 1)**(sympy.S(3)/2)
    F = -sqrt(cot(x)**2 - 1)*cot(x)/2 + 5*atanh(cot(x)/sqrt(cot(x)**2 - 1))/2 - 2*sqrt(2)*atanh(sqrt(2)*cot(x)/sqrt(cot(x)**2 - 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_42():
    f = sqrt(cot(x)**2 - 1)
    F = -atanh(cot(x)/sqrt(cot(x)**2 - 1)) + sqrt(2)*atanh(sqrt(2)*cot(x)/sqrt(cot(x)**2 - 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_43():
    f = 1/sqrt(cot(x)**2 - 1)
    F = -sqrt(2)*atanh(sqrt(2)*cot(x)/sqrt(cot(x)**2 - 1))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_44():
    f = cot(x)**3/sqrt(a + b*cot(x)**2)
    F = -atanh(sqrt(a + b*cot(x)**2)/sqrt(a - b))/sqrt(a - b) - sqrt(a + b*cot(x)**2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_45():
    f = cot(x)**2/sqrt(a + b*cot(x)**2)
    F = atan(sqrt(a - b)*cot(x)/sqrt(a + b*cot(x)**2))/sqrt(a - b) - atanh(sqrt(b)*cot(x)/sqrt(a + b*cot(x)**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_46():
    f = cot(x)/sqrt(a + b*cot(x)**2)
    F = atanh(sqrt(a + b*cot(x)**2)/sqrt(a - b))/sqrt(a - b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_47():
    f = tan(x)/sqrt(a + b*cot(x)**2)
    F = -atanh(sqrt(a + b*cot(x)**2)/sqrt(a - b))/sqrt(a - b) + atanh(sqrt(a + b*cot(x)**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_48():
    f = tan(x)**2/sqrt(a + b*cot(x)**2)
    F = atan(sqrt(a - b)*cot(x)/sqrt(a + b*cot(x)**2))/sqrt(a - b) + sqrt(a + b*cot(x)**2)*tan(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_49():
    f = cot(x)**3/(a + b*cot(x)**2)**(sympy.S(3)/2)
    F = a/(b*(a - b)*sqrt(a + b*cot(x)**2)) - atanh(sqrt(a + b*cot(x)**2)/sqrt(a - b))/(a - b)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_50():
    f = cot(x)**2/(a + b*cot(x)**2)**(sympy.S(3)/2)
    F = -cot(x)/((a - b)*sqrt(a + b*cot(x)**2)) + atan(sqrt(a - b)*cot(x)/sqrt(a + b*cot(x)**2))/(a - b)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_51():
    f = cot(x)/(a + b*cot(x)**2)**(sympy.S(3)/2)
    F = -1/((a - b)*sqrt(a + b*cot(x)**2)) + atanh(sqrt(a + b*cot(x)**2)/sqrt(a - b))/(a - b)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_52():
    f = tan(x)/(a + b*cot(x)**2)**(sympy.S(3)/2)
    F = -atanh(sqrt(a + b*cot(x)**2)/sqrt(a - b))/(a - b)**(sympy.S(3)/2) + b/(a*(a - b)*sqrt(a + b*cot(x)**2)) + atanh(sqrt(a + b*cot(x)**2)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_53():
    f = tan(x)**2/(a + b*cot(x)**2)**(sympy.S(3)/2)
    F = atan(sqrt(a - b)*cot(x)/sqrt(a + b*cot(x)**2))/(a - b)**(sympy.S(3)/2) + b*tan(x)/(a*(a - b)*sqrt(a + b*cot(x)**2)) + (a - 2*b)*sqrt(a + b*cot(x)**2)*tan(x)/(a**2*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_54():
    f = cot(x)**3/(a + b*cot(x)**2)**(sympy.S(5)/2)
    F = a/(b*(a + b*cot(x)**2)**(sympy.S(3)/2)*(3*a - 3*b)) + 1/((a - b)**2*sqrt(a + b*cot(x)**2)) - atanh(sqrt(a + b*cot(x)**2)/sqrt(a - b))/(a - b)**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_55():
    f = cot(x)**2/(a + b*cot(x)**2)**(sympy.S(5)/2)
    F = -cot(x)/((a + b*cot(x)**2)**(sympy.S(3)/2)*(3*a - 3*b)) + atan(sqrt(a - b)*cot(x)/sqrt(a + b*cot(x)**2))/(a - b)**(sympy.S(5)/2) - (2*a + b)*cot(x)/(3*a*(a - b)**2*sqrt(a + b*cot(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_56():
    f = cot(x)/(a + b*cot(x)**2)**(sympy.S(5)/2)
    F = -1/((a + b*cot(x)**2)**(sympy.S(3)/2)*(3*a - 3*b)) - 1/((a - b)**2*sqrt(a + b*cot(x)**2)) + atanh(sqrt(a + b*cot(x)**2)/sqrt(a - b))/(a - b)**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_57():
    f = tan(x)/(a + b*cot(x)**2)**(sympy.S(5)/2)
    F = -atanh(sqrt(a + b*cot(x)**2)/sqrt(a - b))/(a - b)**(sympy.S(5)/2) + b/(3*a*(a - b)*(a + b*cot(x)**2)**(sympy.S(3)/2)) + b*(2*a - b)/(a**2*(a - b)**2*sqrt(a + b*cot(x)**2)) + atanh(sqrt(a + b*cot(x)**2)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_58():
    f = tan(x)**2/(a + b*cot(x)**2)**(sympy.S(5)/2)
    F = atan(sqrt(a - b)*cot(x)/sqrt(a + b*cot(x)**2))/(a - b)**(sympy.S(5)/2) + b*tan(x)/(3*a*(a - b)*(a + b*cot(x)**2)**(sympy.S(3)/2)) + b*(7*a - 4*b)*tan(x)/(3*a**2*(a - b)**2*sqrt(a + b*cot(x)**2)) + (a - 4*b)*sqrt(a + b*cot(x)**2)*(3*a - 2*b)*tan(x)/(3*a**3*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_59():
    f = 1/(cot(x)**3 + 1)
    F = x/2 - log(cot(x) + 1)/6 + log(cot(x)**2 - cot(x) + 1)/3 + log(sin(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_60():
    f = sqrt(a + b*cot(x)**4)*cot(x)
    F = sqrt(b)*atanh(sqrt(b)*cot(x)**2/sqrt(a + b*cot(x)**4))/2 + sqrt(a + b)*atanh((a - b*cot(x)**2)/(sqrt(a + b)*sqrt(a + b*cot(x)**4)))/2 - sqrt(a + b*cot(x)**4)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_61():
    f = (a + b*cot(x)**4)**(sympy.S(3)/2)*cot(x)
    F = sqrt(b)*(3*a + 2*b)*atanh(sqrt(b)*cot(x)**2/sqrt(a + b*cot(x)**4))/4 + (a + b)**(sympy.S(3)/2)*atanh((a - b*cot(x)**2)/(sqrt(a + b)*sqrt(a + b*cot(x)**4)))/2 - (a + b*cot(x)**4)**(sympy.S(3)/2)/6 - sqrt(a + b*cot(x)**4)*(a/2 - b*cot(x)**2/4 + b/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_62():
    f = cot(x)/sqrt(a + b*cot(x)**4)
    F = atanh((a - b*cot(x)**2)/(sqrt(a + b)*sqrt(a + b*cot(x)**4)))/(2*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_63():
    f = cot(x)/(a + b*cot(x)**4)**(sympy.S(3)/2)
    F = atanh((a - b*cot(x)**2)/(sqrt(a + b)*sqrt(a + b*cot(x)**4)))/(2*(a + b)**(sympy.S(3)/2)) - (a + b*cot(x)**2)/(2*a*(a + b)*sqrt(a + b*cot(x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_7_d_trig_pow_m_a_plus_b_c_cot_pow_n_pow_p_64():
    f = cot(x)/(a + b*cot(x)**4)**(sympy.S(5)/2)
    F = atanh((a - b*cot(x)**2)/(sqrt(a + b)*sqrt(a + b*cot(x)**4)))/(2*(a + b)**(sympy.S(5)/2)) - (a + b*cot(x)**2)/(6*a*(a + b)*(a + b*cot(x)**4)**(sympy.S(3)/2)) - (3*a**2 + b*(5*a + 2*b)*cot(x)**2)/(6*a**2*(a + b)**2*sqrt(a + b*cot(x)**4))
    assert integrate(f, x) == F

