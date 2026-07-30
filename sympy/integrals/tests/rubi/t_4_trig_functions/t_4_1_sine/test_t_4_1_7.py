"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.1 Sine/4.1.7 (d trig)^m (a+b (c sin)^n)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n, p = symbols('a b c d e f m n p')

def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_1():
    f = (a*sin(x)**2)**(sympy.S(5)/2)
    F = -8*a**2*sqrt(a*sin(x)**2)*cot(x)/15 - 4*a*(a*sin(x)**2)**(sympy.S(3)/2)*cot(x)/15 - (a*sin(x)**2)**(sympy.S(5)/2)*cot(x)/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_2():
    f = (a*sin(x)**2)**(sympy.S(3)/2)
    F = -2*a*sqrt(a*sin(x)**2)*cot(x)/3 - (a*sin(x)**2)**(sympy.S(3)/2)*cot(x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_3():
    f = sqrt(a*sin(x)**2)
    F = -sqrt(a*sin(x)**2)*cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_4():
    f = 1/sqrt(a*sin(x)**2)
    F = -sin(x)*atanh(cos(x))/sqrt(a*sin(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_5():
    f = (a*sin(x)**2)**(sympy.S(-3)/2)
    F = -sin(x)*atanh(cos(x))/(2*a*sqrt(a*sin(x)**2)) - cot(x)/(2*a*sqrt(a*sin(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_6():
    f = (a*sin(x)**2)**(sympy.S(-5)/2)
    F = -cot(x)/(4*a*(a*sin(x)**2)**(sympy.S(3)/2)) - 3*sin(x)*atanh(cos(x))/(8*a**2*sqrt(a*sin(x)**2)) - 3*cot(x)/(8*a**2*sqrt(a*sin(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_7():
    f = (a*sin(x)**3)**(sympy.S(5)/2)
    F = -2*a**2*sqrt(a*sin(x)**3)*sin(x)**5*cos(x)/15 - 26*a**2*sqrt(a*sin(x)**3)*sin(x)**3*cos(x)/165 - 78*a**2*sqrt(a*sin(x)**3)*sin(x)*cos(x)/385 - 26*a**2*sqrt(a*sin(x)**3)*cot(x)/77 + 26*a**2*sqrt(a*sin(x)**3)*elliptic_f(x/2 - pi/4, 2)/(77*sin(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_8():
    f = (a*sin(x)**3)**(sympy.S(3)/2)
    F = -2*a*sqrt(a*sin(x)**3)*sin(x)**2*cos(x)/9 - 14*a*sqrt(a*sin(x)**3)*cos(x)/45 + 14*a*sqrt(a*sin(x)**3)*elliptic_e(x/2 - pi/4, 2)/(15*sin(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_9():
    f = sqrt(a*sin(x)**3)
    F = -2*sqrt(a*sin(x)**3)*cot(x)/3 + 2*sqrt(a*sin(x)**3)*elliptic_f(x/2 - pi/4, 2)/(3*sin(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_10():
    f = 1/sqrt(a*sin(x)**3)
    F = -2*sin(x)**(sympy.S(3)/2)*elliptic_e(x/2 - pi/4, 2)/sqrt(a*sin(x)**3) - 2*sin(x)*cos(x)/sqrt(a*sin(x)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_11():
    f = (a*sin(x)**3)**(sympy.S(-3)/2)
    F = 10*sin(x)**(sympy.S(3)/2)*elliptic_f(x/2 - pi/4, 2)/(21*a*sqrt(a*sin(x)**3)) - 10*cos(x)/(21*a*sqrt(a*sin(x)**3)) - 2*cot(x)*csc(x)/(7*a*sqrt(a*sin(x)**3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_12():
    f = (a*sin(x)**3)**(sympy.S(-5)/2)
    F = -154*sin(x)**(sympy.S(3)/2)*elliptic_e(x/2 - pi/4, 2)/(195*a**2*sqrt(a*sin(x)**3)) - 154*sin(x)*cos(x)/(195*a**2*sqrt(a*sin(x)**3)) - 2*cot(x)*csc(x)**4/(13*a**2*sqrt(a*sin(x)**3)) - 22*cot(x)*csc(x)**2/(117*a**2*sqrt(a*sin(x)**3)) - 154*cot(x)/(585*a**2*sqrt(a*sin(x)**3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_13():
    f = (a*sin(x)**4)**(sympy.S(5)/2)
    F = 63*a**2*x*sqrt(a*sin(x)**4)*csc(x)**2/256 - a**2*sqrt(a*sin(x)**4)*sin(x)**7*cos(x)/10 - 9*a**2*sqrt(a*sin(x)**4)*sin(x)**5*cos(x)/80 - 21*a**2*sqrt(a*sin(x)**4)*sin(x)**3*cos(x)/160 - 21*a**2*sqrt(a*sin(x)**4)*sin(x)*cos(x)/128 - 63*a**2*sqrt(a*sin(x)**4)*cot(x)/256
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_14():
    f = (a*sin(x)**4)**(sympy.S(3)/2)
    F = 5*a*x*sqrt(a*sin(x)**4)*csc(x)**2/16 - a*sqrt(a*sin(x)**4)*sin(x)**3*cos(x)/6 - 5*a*sqrt(a*sin(x)**4)*sin(x)*cos(x)/24 - 5*a*sqrt(a*sin(x)**4)*cot(x)/16
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_15():
    f = sqrt(a*sin(x)**4)
    F = x*sqrt(a*sin(x)**4)*csc(x)**2/2 - sqrt(a*sin(x)**4)*cot(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_16():
    f = 1/sqrt(a*sin(x)**4)
    F = -sin(x)*cos(x)/sqrt(a*sin(x)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_17():
    f = (a*sin(x)**4)**(sympy.S(-3)/2)
    F = -sin(x)*cos(x)/(a*sqrt(a*sin(x)**4)) - cos(x)**2*cot(x)**3/(5*a*sqrt(a*sin(x)**4)) - 2*cos(x)**2*cot(x)/(3*a*sqrt(a*sin(x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_18():
    f = (a*sin(x)**4)**(sympy.S(-5)/2)
    F = -sin(x)*cos(x)/(a**2*sqrt(a*sin(x)**4)) - cos(x)**2*cot(x)**7/(9*a**2*sqrt(a*sin(x)**4)) - 4*cos(x)**2*cot(x)**5/(7*a**2*sqrt(a*sin(x)**4)) - 6*cos(x)**2*cot(x)**3/(5*a**2*sqrt(a*sin(x)**4)) - 4*cos(x)**2*cot(x)/(3*a**2*sqrt(a*sin(x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_19():
    f = (c*sin(a + b*x)**m)**(sympy.S(5)/2)
    F = 2*c**2*sqrt(c*sin(a + b*x)**m)*sin(a + b*x)**(2*m + 1)*cos(a + b*x)*hyper((sympy.S.Half, 5*m/4 + sympy.S.Half), (5*m/4 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*(5*m + 2)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_20():
    f = (c*sin(a + b*x)**m)**(sympy.S(3)/2)
    F = 2*c*sqrt(c*sin(a + b*x)**m)*sin(a + b*x)**(m + 1)*cos(a + b*x)*hyper((sympy.S.Half, 3*m/4 + sympy.S.Half), (3*m/4 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*(3*m + 2)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_21():
    f = sqrt(c*sin(a + b*x)**m)
    F = 2*sqrt(c*sin(a + b*x)**m)*sin(a + b*x)*cos(a + b*x)*hyper((sympy.S.Half, m/4 + sympy.S.Half), (m/4 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*(m + 2)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_22():
    f = 1/sqrt(c*sin(a + b*x)**m)
    F = 2*sin(a + b*x)*cos(a + b*x)*hyper((sympy.S.Half, sympy.S.Half - m/4), (sympy.S(3)/2 - m/4,), sin(a + b*x)**2)/(b*sqrt(c*sin(a + b*x)**m)*(2 - m)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_23():
    f = (c*sin(a + b*x)**m)**(sympy.S(-3)/2)
    F = 2*sin(a + b*x)**(1 - m)*cos(a + b*x)*hyper((sympy.S.Half, sympy.S.Half - 3*m/4), (sympy.S(3)/2 - 3*m/4,), sin(a + b*x)**2)/(b*c*sqrt(c*sin(a + b*x)**m)*(2 - 3*m)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_24():
    f = (c*sin(a + b*x)**m)**(sympy.S(-5)/2)
    F = 2*sin(a + b*x)**(1 - 2*m)*cos(a + b*x)*hyper((sympy.S.Half, sympy.S.Half - 5*m/4), (sympy.S(3)/2 - 5*m/4,), sin(a + b*x)**2)/(b*c**2*sqrt(c*sin(a + b*x)**m)*(2 - 5*m)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_25():
    f = (b*sin(c + d*x)**n)**p
    F = (b*sin(c + d*x)**n)**p*sin(c + d*x)*cos(c + d*x)*hyper((sympy.S.Half, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*(n*p + 1)*sqrt(cos(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_26():
    f = (c*sin(a + b*x)**2)**p
    F = (c*sin(a + b*x)**2)**p*sin(a + b*x)*cos(a + b*x)*hyper((sympy.S.Half, p + sympy.S.Half), (p + sympy.S(3)/2,), sin(a + b*x)**2)/(b*(2*p + 1)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_27():
    f = (c*sin(a + b*x)**3)**p
    F = (c*sin(a + b*x)**3)**p*sin(a + b*x)*cos(a + b*x)*hyper((sympy.S.Half, 3*p/2 + sympy.S.Half), (3*p/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*(3*p + 1)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_28():
    f = (c*sin(a + b*x)**4)**p
    F = (c*sin(a + b*x)**4)**p*sin(a + b*x)*cos(a + b*x)*hyper((sympy.S.Half, 2*p + sympy.S.Half), (2*p + sympy.S(3)/2,), sin(a + b*x)**2)/(b*(4*p + 1)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_29():
    f = (c*sin(a + b*x)**n)**(1/n)
    F = -(c*sin(a + b*x)**n)**(1/n)*cot(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_30():
    f = (a*(b*sin(c + d*x))**p)**n
    F = (a*(b*sin(c + d*x))**p)**n*sin(c + d*x)*cos(c + d*x)*hyper((sympy.S.Half, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*(n*p + 1)*sqrt(cos(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_31():
    f = -a*sin(x)**2 + a
    F = a*x/2 + a*sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_32():
    f = (-a*sin(x)**2 + a)**2
    F = 3*a**2*x/8 + a**2*sin(x)*cos(x)**3/4 + 3*a**2*sin(x)*cos(x)/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_33():
    f = (-a*sin(x)**2 + a)**3
    F = 5*a**3*x/16 + a**3*sin(x)*cos(x)**5/6 + 5*a**3*sin(x)*cos(x)**3/24 + 5*a**3*sin(x)*cos(x)/16
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_34():
    f = (-a*sin(x)**2 + a)**4
    F = 35*a**4*x/128 + a**4*sin(x)*cos(x)**7/8 + 7*a**4*sin(x)*cos(x)**5/48 + 35*a**4*sin(x)*cos(x)**3/192 + 35*a**4*sin(x)*cos(x)/128
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_35():
    f = sin(c + d*x)**7/(-a*sin(c + d*x)**2 + a)
    F = cos(c + d*x)**5/(5*a*d) - cos(c + d*x)**3/(a*d) + 3*cos(c + d*x)/(a*d) + sec(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_36():
    f = sin(c + d*x)**5/(-a*sin(c + d*x)**2 + a)
    F = -cos(c + d*x)**3/(3*a*d) + 2*cos(c + d*x)/(a*d) + sec(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_37():
    f = sin(c + d*x)**3/(-a*sin(c + d*x)**2 + a)
    F = cos(c + d*x)/(a*d) + sec(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_38():
    f = sin(c + d*x)/(-a*sin(c + d*x)**2 + a)
    F = sec(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_39():
    f = csc(c + d*x)/(-a*sin(c + d*x)**2 + a)
    F = -atanh(cos(c + d*x))/(a*d) + sec(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_40():
    f = csc(c + d*x)**3/(-a*sin(c + d*x)**2 + a)
    F = -3*atanh(cos(c + d*x))/(2*a*d) - csc(c + d*x)**2*sec(c + d*x)/(2*a*d) + 3*sec(c + d*x)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_41():
    f = csc(c + d*x)**5/(-a*sin(c + d*x)**2 + a)
    F = -15*atanh(cos(c + d*x))/(8*a*d) - csc(c + d*x)**4*sec(c + d*x)/(4*a*d) - 5*csc(c + d*x)**2*sec(c + d*x)/(8*a*d) + 15*sec(c + d*x)/(8*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_42():
    f = sin(c + d*x)**6/(-a*sin(c + d*x)**2 + a)
    F = -15*x/(8*a) - sin(c + d*x)**4*tan(c + d*x)/(4*a*d) - 5*sin(c + d*x)**2*tan(c + d*x)/(8*a*d) + 15*tan(c + d*x)/(8*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_43():
    f = sin(c + d*x)**4/(-a*sin(c + d*x)**2 + a)
    F = -3*x/(2*a) - sin(c + d*x)**2*tan(c + d*x)/(2*a*d) + 3*tan(c + d*x)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_44():
    f = sin(c + d*x)**2/(-a*sin(c + d*x)**2 + a)
    F = -x/a + tan(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_45():
    f = 1/(-a*sin(c + d*x)**2 + a)
    F = tan(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_46():
    f = csc(c + d*x)**2/(-a*sin(c + d*x)**2 + a)
    F = tan(c + d*x)/(a*d) - cot(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_47():
    f = csc(c + d*x)**4/(-a*sin(c + d*x)**2 + a)
    F = tan(c + d*x)/(a*d) - cot(c + d*x)**3/(3*a*d) - 2*cot(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_48():
    f = csc(c + d*x)**6/(-a*sin(c + d*x)**2 + a)
    F = tan(c + d*x)/(a*d) - cot(c + d*x)**5/(5*a*d) - cot(c + d*x)**3/(a*d) - 3*cot(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_49():
    f = sin(c + d*x)**7/(-a*sin(c + d*x)**2 + a)**2
    F = cos(c + d*x)**3/(3*a**2*d) - 3*cos(c + d*x)/(a**2*d) + sec(c + d*x)**3/(3*a**2*d) - 3*sec(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_50():
    f = sin(c + d*x)**5/(-a*sin(c + d*x)**2 + a)**2
    F = -cos(c + d*x)/(a**2*d) + sec(c + d*x)**3/(3*a**2*d) - 2*sec(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_51():
    f = sin(c + d*x)**3/(-a*sin(c + d*x)**2 + a)**2
    F = sec(c + d*x)**3/(3*a**2*d) - sec(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_52():
    f = sin(c + d*x)/(-a*sin(c + d*x)**2 + a)**2
    F = sec(c + d*x)**3/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_53():
    f = csc(c + d*x)/(-a*sin(c + d*x)**2 + a)**2
    F = -atanh(cos(c + d*x))/(a**2*d) + sec(c + d*x)**3/(3*a**2*d) + sec(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_54():
    f = csc(c + d*x)**3/(-a*sin(c + d*x)**2 + a)**2
    F = -5*atanh(cos(c + d*x))/(2*a**2*d) - csc(c + d*x)**2*sec(c + d*x)**3/(2*a**2*d) + 5*sec(c + d*x)**3/(6*a**2*d) + 5*sec(c + d*x)/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_55():
    f = sin(c + d*x)**6/(-a*sin(c + d*x)**2 + a)**2
    F = 5*x/(2*a**2) - sin(c + d*x)**2*tan(c + d*x)**3/(2*a**2*d) + 5*tan(c + d*x)**3/(6*a**2*d) - 5*tan(c + d*x)/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_56():
    f = sin(c + d*x)**4/(-a*sin(c + d*x)**2 + a)**2
    F = x/a**2 + tan(c + d*x)**3/(3*a**2*d) - tan(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_57():
    f = sin(c + d*x)**2/(-a*sin(c + d*x)**2 + a)**2
    F = tan(c + d*x)**3/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_58():
    f = (-a*sin(c + d*x)**2 + a)**(-2)
    F = tan(c + d*x)**3/(3*a**2*d) + tan(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_59():
    f = csc(c + d*x)**2/(-a*sin(c + d*x)**2 + a)**2
    F = tan(c + d*x)**3/(3*a**2*d) + 2*tan(c + d*x)/(a**2*d) - cot(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_60():
    f = csc(c + d*x)**4/(-a*sin(c + d*x)**2 + a)**2
    F = tan(c + d*x)**3/(3*a**2*d) + 3*tan(c + d*x)/(a**2*d) - cot(c + d*x)**3/(3*a**2*d) - 3*cot(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_61():
    f = (-a*sin(x)**2 + a)**(-3)
    F = tan(x)**5/(5*a**3) + 2*tan(x)**3/(3*a**3) + tan(x)/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_62():
    f = (-a*sin(x)**2 + a)**(-4)
    F = tan(x)**7/(7*a**4) + 3*tan(x)**5/(5*a**4) + tan(x)**3/a**4 + tan(x)/a**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_63():
    f = (-a*sin(x)**2 + a)**(-5)
    F = tan(x)**9/(9*a**5) + 4*tan(x)**7/(7*a**5) + 6*tan(x)**5/(5*a**5) + 4*tan(x)**3/(3*a**5) + tan(x)/a**5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_64():
    f = (a + b*sin(c + d*x)**2)*sin(c + d*x)**3
    F = -b*cos(c + d*x)**5/(5*d) - (a + b)*cos(c + d*x)/d + (a + 2*b)*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_65():
    f = (a + b*sin(c + d*x)**2)*sin(c + d*x)
    F = b*cos(c + d*x)**3/(3*d) - (a + b)*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_66():
    f = (a + b*sin(c + d*x)**2)*csc(c + d*x)
    F = -a*atanh(cos(c + d*x))/d - b*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_67():
    f = (a + b*sin(c + d*x)**2)*csc(c + d*x)**3
    F = -a*cot(c + d*x)*csc(c + d*x)/(2*d) - (a + 2*b)*atanh(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_68():
    f = (a + b*sin(c + d*x)**2)*sin(c + d*x)**4
    F = -b*sin(c + d*x)**5*cos(c + d*x)/(6*d) + x*(3*a/8 + 5*b/16) - (6*a + 5*b)*sin(c + d*x)**3*cos(c + d*x)/(24*d) - (6*a + 5*b)*sin(c + d*x)*cos(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_69():
    f = (a + b*sin(c + d*x)**2)*sin(c + d*x)**2
    F = -b*sin(c + d*x)**3*cos(c + d*x)/(4*d) + x*(a/2 + 3*b/8) - (4*a + 3*b)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_70():
    f = a + b*sin(c + d*x)**2
    F = a*x + b*x/2 - b*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_71():
    f = (a + b*sin(c + d*x)**2)*csc(c + d*x)**2
    F = -a*cot(c + d*x)/d + b*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_72():
    f = (a + b*sin(c + d*x)**2)*csc(c + d*x)**4
    F = -a*cot(c + d*x)*csc(c + d*x)**2/(3*d) - (2*a + 3*b)*cot(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_73():
    f = (a + b*sin(c + d*x)**2)*csc(c + d*x)**6
    F = -a*cot(c + d*x)*csc(c + d*x)**4/(5*d) - (4*a + 5*b)*cot(c + d*x)**3/(15*d) - (4*a + 5*b)*cot(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_74():
    f = a + b*sin(x)**2
    F = a*x + b*x/2 - b*sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_75():
    f = (a + b*sin(x)**2)**2
    F = -b**2*sin(x)**3*cos(x)/4 - b*(8*a + 3*b)*sin(x)*cos(x)/8 + x*(a**2 + a*b + 3*b**2/8)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_76():
    f = (a + b*sin(x)**2)**3
    F = -5*b**2*(2*a + b)*sin(x)**3*cos(x)/24 - b*(a + b*sin(x)**2)**2*sin(x)*cos(x)/6 - b*(64*a**2 + 54*a*b + 15*b**2)*sin(x)*cos(x)/48 + x*(a/8 + b/16)*(8*a**2 + 8*a*b + 5*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_77():
    f = (a + b*sin(x)**2)**4
    F = -b**2*(104*a**2 + 104*a*b + 35*b**2)*sin(x)**3*cos(x)/192 - b*(a + b*sin(x)**2)**3*sin(x)*cos(x)/8 - 7*b*(a + b*sin(x)**2)**2*(2*a + b)*sin(x)*cos(x)/48 - b*(608*a**3 + 808*a**2*b + 480*a*b**2 + 105*b**3)*sin(x)*cos(x)/384 + x*(a**4 + 2*a**3*b + 9*a**2*b**2/4 + 5*a*b**3/4 + 35*b**4/128)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_78():
    f = sin(c + d*x)**7/(a + b*sin(c + d*x)**2)
    F = a**3*atanh(sqrt(b)*cos(c + d*x)/sqrt(a + b))/(b**(sympy.S(7)/2)*d*sqrt(a + b)) - cos(c + d*x)**5/(5*b*d) - (a - 2*b)*cos(c + d*x)**3/(3*b**2*d) - (a**2 - a*b + b**2)*cos(c + d*x)/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_79():
    f = sin(c + d*x)**5/(a + b*sin(c + d*x)**2)
    F = -a**2*atanh(sqrt(b)*cos(c + d*x)/sqrt(a + b))/(b**(sympy.S(5)/2)*d*sqrt(a + b)) + cos(c + d*x)**3/(3*b*d) + (a - b)*cos(c + d*x)/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_80():
    f = sin(c + d*x)**3/(a + b*sin(c + d*x)**2)
    F = a*atanh(sqrt(b)*cos(c + d*x)/sqrt(a + b))/(b**(sympy.S(3)/2)*d*sqrt(a + b)) - cos(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_81():
    f = sin(c + d*x)/(a + b*sin(c + d*x)**2)
    F = -atanh(sqrt(b)*cos(c + d*x)/sqrt(a + b))/(sqrt(b)*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_82():
    f = csc(c + d*x)/(a + b*sin(c + d*x)**2)
    F = sqrt(b)*atanh(sqrt(b)*cos(c + d*x)/sqrt(a + b))/(a*d*sqrt(a + b)) - atanh(cos(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_83():
    f = csc(c + d*x)**3/(a + b*sin(c + d*x)**2)
    F = -cot(c + d*x)*csc(c + d*x)/(2*a*d) - b**(sympy.S(3)/2)*atanh(sqrt(b)*cos(c + d*x)/sqrt(a + b))/(a**2*d*sqrt(a + b)) - (a - 2*b)*atanh(cos(c + d*x))/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_84():
    f = csc(c + d*x)**5/(a + b*sin(c + d*x)**2)
    F = -cot(c + d*x)*csc(c + d*x)**3/(4*a*d) - (3*a - 4*b)*cot(c + d*x)*csc(c + d*x)/(8*a**2*d) + b**(sympy.S(5)/2)*atanh(sqrt(b)*cos(c + d*x)/sqrt(a + b))/(a**3*d*sqrt(a + b)) - (3*a**2 - 4*a*b + 8*b**2)*atanh(cos(c + d*x))/(8*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_85():
    f = sin(c + d*x)**8/(a + b*sin(c + d*x)**2)
    F = a**(sympy.S(7)/2)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(b**4*d*sqrt(a + b)) - sin(c + d*x)**5*cos(c + d*x)/(6*b*d) + (6*a - 5*b)*sin(c + d*x)**3*cos(c + d*x)/(24*b**2*d) - (8*a**2 - 6*a*b + 5*b**2)*sin(c + d*x)*cos(c + d*x)/(16*b**3*d) - x*(16*a**3 - 8*a**2*b + 6*a*b**2 - 5*b**3)/(16*b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_86():
    f = sin(c + d*x)**6/(a + b*sin(c + d*x)**2)
    F = -a**(sympy.S(5)/2)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(b**3*d*sqrt(a + b)) - sin(c + d*x)**3*cos(c + d*x)/(4*b*d) + (4*a - 3*b)*sin(c + d*x)*cos(c + d*x)/(8*b**2*d) + x*(8*a**2 - 4*a*b + 3*b**2)/(8*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_87():
    f = sin(c + d*x)**4/(a + b*sin(c + d*x)**2)
    F = a**(sympy.S(3)/2)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(b**2*d*sqrt(a + b)) - sin(c + d*x)*cos(c + d*x)/(2*b*d) - x*(2*a - b)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_88():
    f = sin(c + d*x)**2/(a + b*sin(c + d*x)**2)
    F = -sqrt(a)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(b*d*sqrt(a + b)) + x/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_89():
    f = 1/(a + b*sin(c + d*x)**2)
    F = atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(sqrt(a)*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_90():
    f = csc(c + d*x)**2/(a + b*sin(c + d*x)**2)
    F = -cot(c + d*x)/(a*d) - b*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(a**(sympy.S(3)/2)*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_91():
    f = csc(c + d*x)**4/(a + b*sin(c + d*x)**2)
    F = -cot(c + d*x)**3/(3*a*d) - (a - b)*cot(c + d*x)/(a**2*d) + b**2*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(a**(sympy.S(5)/2)*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_92():
    f = csc(c + d*x)**6/(a + b*sin(c + d*x)**2)
    F = -cot(c + d*x)**5/(5*a*d) - (2*a - b)*cot(c + d*x)**3/(3*a**2*d) - (a**2 - a*b + b**2)*cot(c + d*x)/(a**3*d) - b**3*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(a**(sympy.S(7)/2)*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_93():
    f = csc(c + d*x)**8/(a + b*sin(c + d*x)**2)
    F = -cot(c + d*x)**7/(7*a*d) - (3*a - b)*cot(c + d*x)**5/(5*a**2*d) - (3*a**2 - 2*a*b + b**2)*cot(c + d*x)**3/(3*a**3*d) - (a - b)*(a**2 + b**2)*cot(c + d*x)/(a**4*d) + b**4*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(a**(sympy.S(9)/2)*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_94():
    f = sin(c + d*x)**7/(a + b*sin(c + d*x)**2)**2
    F = a**3*cos(c + d*x)/(2*b**3*d*(a + b)*(a - b*cos(c + d*x)**2 + b)) - a**2*(5*a + 6*b)*atanh(sqrt(b)*cos(c + d*x)/sqrt(a + b))/(2*b**(sympy.S(7)/2)*d*(a + b)**(sympy.S(3)/2)) + cos(c + d*x)**3/(3*b**2*d) + (2*a - b)*cos(c + d*x)/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_95():
    f = sin(c + d*x)**5/(a + b*sin(c + d*x)**2)**2
    F = -a**2*cos(c + d*x)/(2*b**2*d*(a + b)*(a - b*cos(c + d*x)**2 + b)) + a*(3*a + 4*b)*atanh(sqrt(b)*cos(c + d*x)/sqrt(a + b))/(2*b**(sympy.S(5)/2)*d*(a + b)**(sympy.S(3)/2)) - cos(c + d*x)/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_96():
    f = sin(c + d*x)**3/(a + b*sin(c + d*x)**2)**2
    F = a*cos(c + d*x)/(2*b*d*(a + b)*(a - b*cos(c + d*x)**2 + b)) - (a + 2*b)*atanh(sqrt(b)*cos(c + d*x)/sqrt(a + b))/(2*b**(sympy.S(3)/2)*d*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_97():
    f = sin(c + d*x)/(a + b*sin(c + d*x)**2)**2
    F = -cos(c + d*x)/(d*(2*a + 2*b)*(a - b*cos(c + d*x)**2 + b)) - atanh(sqrt(b)*cos(c + d*x)/sqrt(a + b))/(2*sqrt(b)*d*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_98():
    f = csc(c + d*x)/(a + b*sin(c + d*x)**2)**2
    F = b*cos(c + d*x)/(2*a*d*(a + b)*(a - b*cos(c + d*x)**2 + b)) + sqrt(b)*(3*a + 2*b)*atanh(sqrt(b)*cos(c + d*x)/sqrt(a + b))/(2*a**2*d*(a + b)**(sympy.S(3)/2)) - atanh(cos(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_99():
    f = csc(c + d*x)**3/(a + b*sin(c + d*x)**2)**2
    F = -cot(c + d*x)*csc(c + d*x)/(2*a*d*(a - b*cos(c + d*x)**2 + b)) - b*(a + 2*b)*cos(c + d*x)/(2*a**2*d*(a + b)*(a - b*cos(c + d*x)**2 + b)) - b**(sympy.S(3)/2)*(5*a + 4*b)*atanh(sqrt(b)*cos(c + d*x)/sqrt(a + b))/(2*a**3*d*(a + b)**(sympy.S(3)/2)) - (a - 4*b)*atanh(cos(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_100():
    f = sin(c + d*x)**6/(a + b*sin(c + d*x)**2)**2
    F = a**(sympy.S(3)/2)*(4*a + 5*b)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(2*b**3*d*(a + b)**(sympy.S(3)/2)) - a*(2*a + b)*tan(c + d*x)/(2*b**2*d*(a + b)*(a + (a + b)*tan(c + d*x)**2)) - sin(c + d*x)**2*tan(c + d*x)/(2*b*d*(a + (a + b)*tan(c + d*x)**2)) - x*(4*a - b)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_101():
    f = sin(c + d*x)**4/(a + b*sin(c + d*x)**2)**2
    F = -sqrt(a)*(2*a + 3*b)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(2*b**2*d*(a + b)**(sympy.S(3)/2)) + a*tan(c + d*x)/(2*b*d*(a + b)*(a + (a + b)*tan(c + d*x)**2)) + x/b**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_102():
    f = sin(c + d*x)**2/(a + b*sin(c + d*x)**2)**2
    F = -sin(c + d*x)*cos(c + d*x)/(d*(a + b*sin(c + d*x)**2)*(2*a + 2*b)) + atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(2*sqrt(a)*d*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_103():
    f = (a + b*sin(c + d*x)**2)**(-2)
    F = b*sin(c + d*x)*cos(c + d*x)/(2*a*d*(a + b)*(a + b*sin(c + d*x)**2)) + (2*a + b)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_104():
    f = csc(c + d*x)**4/(a + b*sin(c + d*x)**2)**2
    F = b*csc(c + d*x)**3*sec(c + d*x)/(2*a*d*(a + b)*(a + (a + b)*tan(c + d*x)**2)) - (2*a + 5*b)*cot(c + d*x)**3/(6*a**2*d*(a + b)) - (2*a**2 - a*b - 5*b**2)*cot(c + d*x)/(2*a**3*d*(a + b)) + b**2*(6*a + 5*b)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(2*a**(sympy.S(7)/2)*d*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_105():
    f = sin(c + d*x)**6/(a + b*sin(c + d*x)**2)**3
    F = -sqrt(a)*(8*a**2 + 20*a*b + 15*b**2)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(8*b**3*d*(a + b)**(sympy.S(5)/2)) + a*tan(c + d*x)**3/(4*b*d*(a + b)*(a + (a + b)*tan(c + d*x)**2)**2) + a*(4*a + 7*b)*tan(c + d*x)/(8*b**2*d*(a + b)**2*(a + (a + b)*tan(c + d*x)**2)) + x/b**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_106():
    f = sin(c + d*x)**4/(a + b*sin(c + d*x)**2)**3
    F = -tan(c + d*x)**3/(d*(a + (a + b)*tan(c + d*x)**2)**2*(4*a + 4*b)) - 3*tan(c + d*x)/(8*d*(a + b)**2*(a + (a + b)*tan(c + d*x)**2)) + 3*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(8*sqrt(a)*d*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_107():
    f = sin(c + d*x)**2/(a + b*sin(c + d*x)**2)**3
    F = -sin(c + d*x)*cos(c + d*x)/(d*(a + b*sin(c + d*x)**2)**2*(4*a + 4*b)) - (2*a - b)*sin(c + d*x)*cos(c + d*x)/(8*a*d*(a + b)**2*(a + b*sin(c + d*x)**2)) + (4*a + b)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(8*a**(sympy.S(3)/2)*d*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_108():
    f = (a + b*sin(c + d*x)**2)**(-3)
    F = b*sin(c + d*x)*cos(c + d*x)/(4*a*d*(a + b)*(a + b*sin(c + d*x)**2)**2) + 3*b*(2*a + b)*sin(c + d*x)*cos(c + d*x)/(8*a**2*d*(a + b)**2*(a + b*sin(c + d*x)**2)) + (8*a**2 + 8*a*b + 3*b**2)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_109():
    f = csc(c + d*x)**2/(a + b*sin(c + d*x)**2)**3
    F = b*csc(c + d*x)*sec(c + d*x)**3/(4*a*d*(a + b)*(a + (a + b)*tan(c + d*x)**2)**2) + b*(4*a + 5*b + (4*a + b)*tan(c + d*x)**2)*cot(c + d*x)/(8*a**2*d*(a + b)**2*(a + (a + b)*tan(c + d*x)**2)) - (2*a + 3*b)*(4*a + 5*b)*cot(c + d*x)/(8*a**3*d*(a + b)**2) - 3*b*(8*a**2 + 12*a*b + 5*b**2)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(8*a**(sympy.S(7)/2)*d*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_110():
    f = (a + b*sin(c + d*x)**2)**(-4)
    F = b*sin(c + d*x)*cos(c + d*x)/(6*a*d*(a + b)*(a + b*sin(c + d*x)**2)**3) + 5*b*(2*a + b)*sin(c + d*x)*cos(c + d*x)/(24*a**2*d*(a + b)**2*(a + b*sin(c + d*x)**2)**2) + b*(44*a**2 + 44*a*b + 15*b**2)*sin(c + d*x)*cos(c + d*x)/(48*a**3*d*(a + b)**3*(a + b*sin(c + d*x)**2)) + (2*a + b)*(8*a**2 + 8*a*b + 5*b**2)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(16*a**(sympy.S(7)/2)*d*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_111():
    f = (a + b*sin(c + d*x)**2)**(-5)
    F = b*sin(c + d*x)*cos(c + d*x)/(8*a*d*(a + b)*(a + b*sin(c + d*x)**2)**4) + 7*b*(2*a + b)*sin(c + d*x)*cos(c + d*x)/(48*a**2*d*(a + b)**2*(a + b*sin(c + d*x)**2)**3) + b*(104*a**2 + 104*a*b + 35*b**2)*sin(c + d*x)*cos(c + d*x)/(192*a**3*d*(a + b)**3*(a + b*sin(c + d*x)**2)**2) + 5*b*(2*a + b)*(40*a**2 + 40*a*b + 21*b**2)*sin(c + d*x)*cos(c + d*x)/(384*a**4*d*(a + b)**4*(a + b*sin(c + d*x)**2)) + (128*a**4 + 256*a**3*b + 288*a**2*b**2 + 160*a*b**3 + 35*b**4)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(128*a**(sympy.S(9)/2)*d*(a + b)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_112():
    f = sin(x)/sqrt(sin(x)**2 + 1)
    F = -asin(sqrt(2)*cos(x)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_113():
    f = sqrt(sin(x)**2 + 1)*sin(x)
    F = -sqrt(2 - cos(x)**2)*cos(x)/2 - asin(sqrt(2)*cos(x)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_114():
    f = sin(3*x + 7)/sqrt(sin(3*x + 7)**2 + 3)
    F = -asin(cos(3*x + 7)/2)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_115():
    f = (-a*sin(x)**2 + a)**(sympy.S(5)/2)
    F = 8*a**2*sqrt(a*cos(x)**2)*tan(x)/15 + 4*a*(a*cos(x)**2)**(sympy.S(3)/2)*tan(x)/15 + (a*cos(x)**2)**(sympy.S(5)/2)*tan(x)/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_116():
    f = (-a*sin(x)**2 + a)**(sympy.S(3)/2)
    F = 2*a*sqrt(a*cos(x)**2)*tan(x)/3 + (a*cos(x)**2)**(sympy.S(3)/2)*tan(x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_117():
    f = sqrt(-a*sin(x)**2 + a)
    F = sqrt(a*cos(x)**2)*tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_118():
    f = 1/sqrt(-a*sin(x)**2 + a)
    F = cos(x)*atanh(sin(x))/sqrt(a*cos(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_119():
    f = (-a*sin(x)**2 + a)**(sympy.S(-3)/2)
    F = cos(x)*atanh(sin(x))/(2*a*sqrt(a*cos(x)**2)) + tan(x)/(2*a*sqrt(a*cos(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_120():
    f = (-a*sin(x)**2 + a)**(sympy.S(-5)/2)
    F = tan(x)/(4*a*(a*cos(x)**2)**(sympy.S(3)/2)) + 3*cos(x)*atanh(sin(x))/(8*a**2*sqrt(a*cos(x)**2)) + 3*tan(x)/(8*a**2*sqrt(a*cos(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_121():
    f = sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)**3
    F = (a - 3*b)*sqrt(a - b*cos(e + f*x)**2 + b)*cos(e + f*x)/(8*b*f) - (a - b*cos(e + f*x)**2 + b)**(sympy.S(3)/2)*cos(e + f*x)/(4*b*f) + (a - 3*b)*(a + b)*atan(sqrt(b)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(8*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_122():
    f = sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)
    F = -sqrt(a - b*cos(e + f*x)**2 + b)*cos(e + f*x)/(2*f) - (a + b)*atan(sqrt(b)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(2*sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_123():
    f = sqrt(a + b*sin(e + f*x)**2)*csc(e + f*x)
    F = -sqrt(a)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/f - sqrt(b)*atan(sqrt(b)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_124():
    f = sqrt(a + b*sin(e + f*x)**2)*csc(e + f*x)**3
    F = -sqrt(a - b*cos(e + f*x)**2 + b)*cot(e + f*x)*csc(e + f*x)/(2*f) - (a + b)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(2*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_125():
    f = sqrt(a + b*sin(e + f*x)**2)*csc(e + f*x)**5
    F = -(3*a - b)*sqrt(a - b*cos(e + f*x)**2 + b)*cot(e + f*x)*csc(e + f*x)/(8*a*f) - (a - b*cos(e + f*x)**2 + b)**(sympy.S(3)/2)*cot(e + f*x)*csc(e + f*x)**3/(4*a*f) - (a + b)*(3*a - b)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(8*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_126():
    f = sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)**4
    F = 2*a*sqrt(1 + b*sin(e + f*x)**2/a)*(a - 2*b)*(a + b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(15*b**2*f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)**3*cos(e + f*x)/(5*f) - (a + 4*b)*sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)*cos(e + f*x)/(15*b*f) - sqrt(a + b*sin(e + f*x)**2)*(2*a**2 - 3*a*b - 8*b**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(15*b**2*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_127():
    f = sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)**2
    F = -a*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*elliptic_f(e + f*x, -b/a)/(3*b*f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)*cos(e + f*x)/(3*f) + (a + 2*b)*sqrt(a + b*sin(e + f*x)**2)*elliptic_e(e + f*x, -b/a)/(3*b*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_128():
    f = sqrt(a + b*sin(e + f*x)**2)
    F = sqrt(a + b*sin(e + f*x)**2)*elliptic_e(e + f*x, -b/a)/(f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_129():
    f = sqrt(a + b*sin(e + f*x)**2)*csc(e + f*x)**2
    F = sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)/f - sqrt(a + b*sin(e + f*x)**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_130():
    f = sqrt(a + b*sin(e + f*x)**2)*csc(e + f*x)**4
    F = sqrt(1 + b*sin(e + f*x)**2/a)*(2*a + 2*b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)*csc(e + f*x)**2/(3*f) - sqrt(a + b*sin(e + f*x)**2)*(2*a + b)*cot(e + f*x)/(3*a*f) - sqrt(a + b*sin(e + f*x)**2)*(2*a + b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_131():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)**3
    F = (a - 5*b)*(a + b)*sqrt(a - b*cos(e + f*x)**2 + b)*cos(e + f*x)/(16*b*f) + (a - 5*b)*(a - b*cos(e + f*x)**2 + b)**(sympy.S(3)/2)*cos(e + f*x)/(24*b*f) - (a - b*cos(e + f*x)**2 + b)**(sympy.S(5)/2)*cos(e + f*x)/(6*b*f) + (a - 5*b)*(a + b)**2*atan(sqrt(b)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(16*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_132():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)
    F = -(3*a + 3*b)*sqrt(a - b*cos(e + f*x)**2 + b)*cos(e + f*x)/(8*f) - (a - b*cos(e + f*x)**2 + b)**(sympy.S(3)/2)*cos(e + f*x)/(4*f) - 3*(a + b)**2*atan(sqrt(b)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(8*sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_133():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/f - sqrt(b)*(3*a + b)*atan(sqrt(b)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(2*f) - b*sqrt(a - b*cos(e + f*x)**2 + b)*cos(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_134():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**3
    F = -sqrt(a)*(a + 3*b)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(2*f) - a*sqrt(a - b*cos(e + f*x)**2 + b)*cot(e + f*x)*csc(e + f*x)/(2*f) - b**(sympy.S(3)/2)*atan(sqrt(b)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_135():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**5
    F = -(3*a + 3*b)*sqrt(a - b*cos(e + f*x)**2 + b)*cot(e + f*x)*csc(e + f*x)/(8*f) - (a - b*cos(e + f*x)**2 + b)**(sympy.S(3)/2)*cot(e + f*x)*csc(e + f*x)**3/(4*f) - 3*(a + b)**2*atanh(sqrt(a)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(8*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_136():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**7
    F = -(a + b)*(5*a - b)*sqrt(a - b*cos(e + f*x)**2 + b)*cot(e + f*x)*csc(e + f*x)/(16*a*f) - (5*a - b)*(a - b*cos(e + f*x)**2 + b)**(sympy.S(3)/2)*cot(e + f*x)*csc(e + f*x)**3/(24*a*f) - (a - b*cos(e + f*x)**2 + b)**(sympy.S(5)/2)*cot(e + f*x)*csc(e + f*x)**5/(6*a*f) - (a + b)**2*(5*a - b)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(16*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_137():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)**4
    F = a*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*(2*a**2 - 5*a*b - 8*b**2)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(35*b**2*f*sqrt(a + b*sin(e + f*x)**2)) - b*sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)**5*cos(e + f*x)/(7*f) - sqrt(a + b*sin(e + f*x)**2)*(8*a + 6*b)*sin(e + f*x)**3*cos(e + f*x)/(35*f) - sqrt(a + b*sin(e + f*x)**2)*(a**2 + 11*a*b + 8*b**2)*sin(e + f*x)*cos(e + f*x)/(35*b*f) - sqrt(a + b*sin(e + f*x)**2)*(2*a + 4*b)*(a**2 - 4*a*b - 4*b**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(35*b**2*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_138():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)**2
    F = -a*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*(3*a + 4*b)*elliptic_f(e + f*x, -b/a)/(15*b*f*sqrt(a + b*sin(e + f*x)**2)) - (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)*cos(e + f*x)/(5*f) - sqrt(a + b*sin(e + f*x)**2)*(3*a + 4*b)*sin(e + f*x)*cos(e + f*x)/(15*f) + sqrt(a + b*sin(e + f*x)**2)*(3*a**2 + 13*a*b + 8*b**2)*elliptic_e(e + f*x, -b/a)/(15*b*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_139():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = -a*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*elliptic_f(e + f*x, -b/a)/(3*f*sqrt(a + b*sin(e + f*x)**2)) - b*sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)*cos(e + f*x)/(3*f) + sqrt(a + b*sin(e + f*x)**2)*(4*a + 2*b)*elliptic_e(e + f*x, -b/a)/(3*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_140():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**2
    F = a*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*sqrt(a + b*sin(e + f*x)**2)) - a*sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)/f - (a - b)*sqrt(a + b*sin(e + f*x)**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_141():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**4
    F = -a*sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)*csc(e + f*x)**2/(3*f) + sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*(2*a + 3*b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(a + b*sin(e + f*x)**2)) + (-2*a - 4*b)*sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)/(3*f) - sqrt(a + b*sin(e + f*x)**2)*(2*a + 4*b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_142():
    f = (a + b*sin(c + d*x)**2)**(sympy.S(5)/2)
    F = -4*a*sqrt(1 + b*sin(c + d*x)**2/a)*(a + b)*(2*a + b)*elliptic_f(c + d*x, -b/a)/(15*d*sqrt(a + b*sin(c + d*x)**2)) - b*(a + b*sin(c + d*x)**2)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)/(5*d) - 4*b*sqrt(a + b*sin(c + d*x)**2)*(2*a + b)*sin(c + d*x)*cos(c + d*x)/(15*d) + sqrt(a + b*sin(c + d*x)**2)*(23*a**2 + 23*a*b + 8*b**2)*elliptic_e(c + d*x, -b/a)/(15*d*sqrt(1 + b*sin(c + d*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_143():
    f = sin(e + f*x)**3/sqrt(a + b*sin(e + f*x)**2)
    F = -sqrt(a - b*cos(e + f*x)**2 + b)*cos(e + f*x)/(2*b*f) + (a - b)*atan(sqrt(b)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(2*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_144():
    f = sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2)
    F = -atan(sqrt(b)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_145():
    f = csc(e + f*x)/sqrt(a + b*sin(e + f*x)**2)
    F = -atanh(sqrt(a)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_146():
    f = csc(e + f*x)**3/sqrt(a + b*sin(e + f*x)**2)
    F = -sqrt(a - b*cos(e + f*x)**2 + b)*cot(e + f*x)*csc(e + f*x)/(2*a*f) - (a - b)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(2*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_147():
    f = sin(e + f*x)**4/sqrt(a + b*sin(e + f*x)**2)
    F = a*sqrt(1 + b*sin(e + f*x)**2/a)*(2*a - b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*b**2*f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)*cos(e + f*x)/(3*b*f) - sqrt(a + b*sin(e + f*x)**2)*(2*a - 2*b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*b**2*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_148():
    f = sin(e + f*x)**2/sqrt(a + b*sin(e + f*x)**2)
    F = -a*sqrt(1 + b*sin(e + f*x)**2/a)*elliptic_f(e + f*x, -b/a)/(b*f*sqrt(a + b*sin(e + f*x)**2)) + sqrt(a + b*sin(e + f*x)**2)*elliptic_e(e + f*x, -b/a)/(b*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_149():
    f = 1/sqrt(a + b*sin(e + f*x)**2)
    F = sqrt(1 + b*sin(e + f*x)**2/a)*elliptic_f(e + f*x, -b/a)/(f*sqrt(a + b*sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_150():
    f = csc(e + f*x)**2/sqrt(a + b*sin(e + f*x)**2)
    F = sqrt(1 + b*sin(e + f*x)**2/a)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)/(a*f) - sqrt(a + b*sin(e + f*x)**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(a*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_151():
    f = csc(e + f*x)**4/sqrt(a + b*sin(e + f*x)**2)
    F = sqrt(1 + b*sin(e + f*x)**2/a)*(2*a - b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a*f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)*csc(e + f*x)**2/(3*a*f) - sqrt(a + b*sin(e + f*x)**2)*(2*a - 2*b)*cot(e + f*x)/(3*a**2*f) - sqrt(a + b*sin(e + f*x)**2)*(2*a - 2*b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a**2*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_152():
    f = sin(e + f*x)**3/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = a*cos(e + f*x)/(b*f*(a + b)*sqrt(a - b*cos(e + f*x)**2 + b)) - atan(sqrt(b)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_153():
    f = sin(e + f*x)/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = -cos(e + f*x)/(f*(a + b)*sqrt(a - b*cos(e + f*x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_154():
    f = csc(e + f*x)/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = b*cos(e + f*x)/(a*f*(a + b)*sqrt(a - b*cos(e + f*x)**2 + b)) - atanh(sqrt(a)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_155():
    f = csc(e + f*x)**3/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = -cot(e + f*x)*csc(e + f*x)/(2*a*f*sqrt(a - b*cos(e + f*x)**2 + b)) - b*(a + 3*b)*cos(e + f*x)/(2*a**2*f*(a + b)*sqrt(a - b*cos(e + f*x)**2 + b)) - (a - 3*b)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(2*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_156():
    f = sin(e + f*x)**6/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = a*sin(e + f*x)**3*cos(e + f*x)/(b*f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) + a*sqrt(1 + b*sin(e + f*x)**2/a)*(8*a - b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*b**3*f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*(4*a + b)*sin(e + f*x)*cos(e + f*x)/(3*b**2*f*(a + b)) - sqrt(a + b*sin(e + f*x)**2)*(8*a**2 + 3*a*b - 2*b**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*b**3*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_157():
    f = sin(e + f*x)**4/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = a*sin(e + f*x)*cos(e + f*x)/(b*f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) - 2*a*sqrt(1 + b*sin(e + f*x)**2/a)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(b**2*f*sqrt(a + b*sin(e + f*x)**2)) + sqrt(a + b*sin(e + f*x)**2)*(2*a + b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(b**2*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_158():
    f = sin(e + f*x)**2/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = -sin(e + f*x)*cos(e + f*x)/(f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) + sqrt(1 + b*sin(e + f*x)**2/a)*elliptic_f(e + f*x, -b/a)/(b*f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*elliptic_e(e + f*x, -b/a)/(b*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_159():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(-3)/2)
    F = b*sin(e + f*x)*cos(e + f*x)/(a*f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) + sqrt(a + b*sin(e + f*x)**2)*elliptic_e(e + f*x, -b/a)/(a*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_160():
    f = csc(e + f*x)**2/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = b*cot(e + f*x)/(a*f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) + sqrt(1 + b*sin(e + f*x)**2/a)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(a*f*sqrt(a + b*sin(e + f*x)**2)) - (a + 2*b)*sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)/(a**2*f*(a + b)) - (a + 2*b)*sqrt(a + b*sin(e + f*x)**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(a**2*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_161():
    f = sin(e + f*x)**5/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = a*sin(e + f*x)**2*cos(e + f*x)/(3*b*f*(a + b)*(a - b*cos(e + f*x)**2 + b)**(sympy.S(3)/2)) + a*(3*a + 5*b)*cos(e + f*x)/(3*b**2*f*(a + b)**2*sqrt(a - b*cos(e + f*x)**2 + b)) - atan(sqrt(b)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_162():
    f = sin(e + f*x)**3/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = -sin(e + f*x)**2*cos(e + f*x)/(f*(3*a + 3*b)*(a - b*cos(e + f*x)**2 + b)**(sympy.S(3)/2)) - 2*cos(e + f*x)/(3*f*(a + b)**2*sqrt(a - b*cos(e + f*x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_163():
    f = sin(e + f*x)/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = -cos(e + f*x)/(f*(3*a + 3*b)*(a - b*cos(e + f*x)**2 + b)**(sympy.S(3)/2)) - 2*cos(e + f*x)/(3*f*(a + b)**2*sqrt(a - b*cos(e + f*x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_164():
    f = csc(e + f*x)/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = b*cos(e + f*x)/(3*a*f*(a + b)*(a - b*cos(e + f*x)**2 + b)**(sympy.S(3)/2)) + b*(5*a + 3*b)*cos(e + f*x)/(3*a**2*f*(a + b)**2*sqrt(a - b*cos(e + f*x)**2 + b)) - atanh(sqrt(a)*cos(e + f*x)/sqrt(a - b*cos(e + f*x)**2 + b))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_165():
    f = sin(e + f*x)**6/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = a*sin(e + f*x)**3*cos(e + f*x)/(3*b*f*(a + b)*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) + 2*a*(2*a + 3*b)*sin(e + f*x)*cos(e + f*x)/(3*b**2*f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2)) - a*sqrt(1 + b*sin(e + f*x)**2/a)*(8*a + 9*b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*b**3*f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) + sqrt(a + b*sin(e + f*x)**2)*(8*a**2 + 13*a*b + 3*b**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*b**3*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_166():
    f = sin(e + f*x)**4/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = a*sin(e + f*x)*cos(e + f*x)/(3*b*f*(a + b)*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) - (2*a + 4*b)*sin(e + f*x)*cos(e + f*x)/(3*b*f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2)) + sqrt(1 + b*sin(e + f*x)**2/a)*(2*a + 3*b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*b**2*f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*(2*a + 4*b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*b**2*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_167():
    f = sin(e + f*x)**2/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = -sin(e + f*x)*cos(e + f*x)/(f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*(3*a + 3*b)) + sqrt(1 + b*sin(e + f*x)**2/a)*elliptic_f(e + f*x, -b/a)/(3*b*f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) - (a - b)*sin(e + f*x)*cos(e + f*x)/(3*a*f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2)) - (a - b)*sqrt(a + b*sin(e + f*x)**2)*elliptic_e(e + f*x, -b/a)/(3*a*b*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_168():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(-5)/2)
    F = b*sin(e + f*x)*cos(e + f*x)/(3*a*f*(a + b)*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) - sqrt(1 + b*sin(e + f*x)**2/a)*elliptic_f(e + f*x, -b/a)/(3*a*f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) + 2*b*(2*a + b)*sin(e + f*x)*cos(e + f*x)/(3*a**2*f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2)) + sqrt(a + b*sin(e + f*x)**2)*(4*a + 2*b)*elliptic_e(e + f*x, -b/a)/(3*a**2*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_169():
    f = csc(e + f*x)**2/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = b*cot(e + f*x)/(3*a*f*(a + b)*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) + 2*b*(3*a + 2*b)*cot(e + f*x)/(3*a**2*f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2)) + sqrt(1 + b*sin(e + f*x)**2/a)*(3*a + 4*b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a**2*f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*(3*a**2 + 13*a*b + 8*b**2)*cot(e + f*x)/(3*a**3*f*(a + b)**2) - sqrt(a + b*sin(e + f*x)**2)*(3*a**2 + 13*a*b + 8*b**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a**3*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_170():
    f = (d*sin(e + f*x))**m*(a + b*sin(e + f*x)**2)**p
    F = -d*(d*sin(e + f*x))**(m - 1)*(a - b*cos(e + f*x)**2 + b)**p*(sin(e + f*x)**2)**(sympy.S.Half - m/2)*cos(e + f*x)*appellf1(sympy.S.Half, -p, sympy.S.Half - m/2, sympy.S(3)/2, b*cos(e + f*x)**2/(a + b), cos(e + f*x)**2)/(f*(-b*cos(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_171():
    f = (a + b*sin(e + f*x)**2)**p*sin(e + f*x)**5
    F = -(a - b*cos(e + f*x)**2 + b)**(p + 1)*sin(e + f*x)**2*cos(e + f*x)/(b*f*(2*p + 5)) + (3*a - 2*b*(p + 2))*(a - b*cos(e + f*x)**2 + b)**(p + 1)*cos(e + f*x)/(b**2*f*(2*p + 3)*(2*p + 5)) - (a - b*cos(e + f*x)**2 + b)**p*(3*a**2 - 4*a*b*(p + 1) + 4*b**2*(p**2 + 3*p + 2))*cos(e + f*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), b*cos(e + f*x)**2/(a + b))/(b**2*f*(2*p + 3)*(2*p + 5)*(-b*cos(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_172():
    f = (a + b*sin(e + f*x)**2)**p*sin(e + f*x)**3
    F = (a - 2*b*(p + 1))*(a - b*cos(e + f*x)**2 + b)**p*cos(e + f*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), b*cos(e + f*x)**2/(a + b))/(b*f*(2*p + 3)*(-b*cos(e + f*x)**2/(a + b) + 1)**p) - (a - b*cos(e + f*x)**2 + b)**(p + 1)*cos(e + f*x)/(b*f*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_173():
    f = (a + b*sin(e + f*x)**2)**p*sin(e + f*x)
    F = -(a - b*cos(e + f*x)**2 + b)**p*cos(e + f*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), b*cos(e + f*x)**2/(a + b))/(f*(-b*cos(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_174():
    f = (a + b*sin(e + f*x)**2)**p*csc(e + f*x)
    F = -(a - b*cos(e + f*x)**2 + b)**p*cos(e + f*x)*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, cos(e + f*x)**2, b*cos(e + f*x)**2/(a + b))/(f*(-b*cos(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_175():
    f = (a + b*sin(e + f*x)**2)**p*csc(e + f*x)**3
    F = -(a - b*cos(e + f*x)**2 + b)**p*cos(e + f*x)*appellf1(sympy.S.Half, 2, -p, sympy.S(3)/2, cos(e + f*x)**2, b*cos(e + f*x)**2/(a + b))/(f*(-b*cos(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_176():
    f = (a + b*sin(e + f*x)**2)**p*csc(e + f*x)**5
    F = -(a - b*cos(e + f*x)**2 + b)**p*cos(e + f*x)*appellf1(sympy.S.Half, 3, -p, sympy.S(3)/2, cos(e + f*x)**2, b*cos(e + f*x)**2/(a + b))/(f*(-b*cos(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_177():
    f = (a + b*sin(e + f*x)**2)**p*sin(e + f*x)**4
    F = (a + b*sin(e + f*x)**2)**p*sqrt(cos(e + f*x)**2)*sin(e + f*x)**4*tan(e + f*x)*appellf1(sympy.S(5)/2, sympy.S.Half, -p, sympy.S(7)/2, sin(e + f*x)**2, -b*sin(e + f*x)**2/a)/(5*f*(1 + b*sin(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_178():
    f = (a + b*sin(e + f*x)**2)**p*sin(e + f*x)**2
    F = (a + b*sin(e + f*x)**2)**p*(sec(e + f*x)**2)**p*tan(e + f*x)**3*appellf1(sympy.S(3)/2, -p, p + 2, sympy.S(5)/2, -(a + b)*tan(e + f*x)**2/a, -tan(e + f*x)**2)/(3*f*(1 + (a + b)*tan(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_179():
    f = (a + b*sin(e + f*x)**2)**p*csc(e + f*x)**2
    F = -(a + b*sin(e + f*x)**2)**p*sqrt(cos(e + f*x)**2)*appellf1(sympy.S(-1)/2, sympy.S.Half, -p, sympy.S.Half, sin(e + f*x)**2, -b*sin(e + f*x)**2/a)*csc(e + f*x)*sec(e + f*x)/(f*(1 + b*sin(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_180():
    f = (a + b*sin(e + f*x)**2)**p*csc(e + f*x)**4
    F = -(a + b*sin(e + f*x)**2)**p*sqrt(cos(e + f*x)**2)*appellf1(sympy.S(-3)/2, sympy.S.Half, -p, sympy.S(-1)/2, sin(e + f*x)**2, -b*sin(e + f*x)**2/a)*csc(e + f*x)**3*sec(e + f*x)/(3*f*(1 + b*sin(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_181():
    f = sin(c + d*x)**7/(a + b*sin(c + d*x)**3)
    F = 2*(-1)**(sympy.S(2)/3)*a**(sympy.S(5)/3)*atan((-a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*b**(sympy.S(7)/3)*d*sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) + 2*(-1)**(sympy.S(1)/3)*a**(sympy.S(5)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(3*b**(sympy.S(7)/3)*d*sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))) - 2*a**(sympy.S(5)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b**(sympy.S(7)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) + a*cos(c + d*x)/(b**2*d) + 3*x/(8*b) - sin(c + d*x)**3*cos(c + d*x)/(4*b*d) - 3*sin(c + d*x)*cos(c + d*x)/(8*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_182():
    f = sin(c + d*x)**5/(a + b*sin(c + d*x)**3)
    F = 2*a*atanh((-(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*b**(sympy.S(5)/3)*d*sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) + 2*a*atanh(((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*b**(sympy.S(5)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) - 2*a*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b**(sympy.S(5)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) + x/(2*b) - sin(c + d*x)*cos(c + d*x)/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_183():
    f = sin(c + d*x)**3/(a + b*sin(c + d*x)**3)
    F = 2*a**(sympy.S(1)/3)*atan((-1)**(sympy.S(1)/3)*((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*b*d*sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) - 2*a**(sympy.S(1)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(3*b*d*sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))) - 2*a**(sympy.S(1)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) + x/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_184():
    f = sin(c + d*x)/(a + b*sin(c + d*x)**3)
    F = 2*(-1)**(sympy.S(2)/3)*atan((-a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*d*sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) + 2*(-1)**(sympy.S(1)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*d*sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))) - 2*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_185():
    f = csc(c + d*x)/(a + b*sin(c + d*x)**3)
    F = 2*b**(sympy.S(1)/3)*atanh((-(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*a*d*sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) + 2*b**(sympy.S(1)/3)*atanh(((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*a*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) - 2*b**(sympy.S(1)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) - atanh(cos(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_186():
    f = csc(c + d*x)**3/(a + b*sin(c + d*x)**3)
    F = -cot(c + d*x)*csc(c + d*x)/(2*a*d) - atanh(cos(c + d*x))/(2*a*d) + 2*b*atan((-1)**(sympy.S(1)/3)*((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(5)/3)*d*sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) - 2*b*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(5)/3)*d*sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))) - 2*b*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**(sympy.S(5)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_187():
    f = csc(c + d*x)**5/(a + b*sin(c + d*x)**3)
    F = -cot(c + d*x)*csc(c + d*x)**3/(4*a*d) - 3*cot(c + d*x)*csc(c + d*x)/(8*a*d) - 3*atanh(cos(c + d*x))/(8*a*d) + b*cot(c + d*x)/(a**2*d) + 2*(-1)**(sympy.S(2)/3)*b**(sympy.S(5)/3)*atan((-a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(7)/3)*d*sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) + 2*(-1)**(sympy.S(1)/3)*b**(sympy.S(5)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(7)/3)*d*sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))) - 2*b**(sympy.S(5)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**(sympy.S(7)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_188():
    f = sin(c + d*x)**6/(a + b*sin(c + d*x)**3)
    F = -2*a**(sympy.S(4)/3)*atan((-1)**(sympy.S(1)/3)*((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*b**2*d*sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) + 2*a**(sympy.S(4)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(3*b**2*d*sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))) + 2*a**(sympy.S(4)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b**2*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) - a*x/b**2 + cos(c + d*x)**3/(3*b*d) - cos(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_189():
    f = sin(c + d*x)**4/(a + b*sin(c + d*x)**3)
    F = -2*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*atan((-a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*b**(sympy.S(4)/3)*d*sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) - 2*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(3*b**(sympy.S(4)/3)*d*sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))) + 2*a**(sympy.S(2)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b**(sympy.S(4)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) - cos(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_190():
    f = sin(c + d*x)**2/(a + b*sin(c + d*x)**3)
    F = -2*atanh((-(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*b**(sympy.S(2)/3)*d*sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) - 2*atanh(((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*b**(sympy.S(2)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) + 2*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_191():
    f = 1/(a + b*sin(c + d*x)**3)
    F = -2*atan((-1)**(sympy.S(1)/3)*((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) + 2*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))) + 2*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_192():
    f = csc(c + d*x)**2/(a + b*sin(c + d*x)**3)
    F = -cot(c + d*x)/(a*d) - 2*(-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)*atan((-a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(4)/3)*d*sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) - 2*(-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(4)/3)*d*sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))) + 2*b**(sympy.S(2)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**(sympy.S(4)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_193():
    f = csc(c + d*x)**4/(a + b*sin(c + d*x)**3)
    F = -cot(c + d*x)**3/(3*a*d) - cot(c + d*x)/(a*d) - 2*b**(sympy.S(4)/3)*atanh((-(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*a**2*d*sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) - 2*b**(sympy.S(4)/3)*atanh(((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*a**2*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) + 2*b**(sympy.S(4)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**2*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) + b*atanh(cos(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_194():
    f = sin(c + d*x)**9/(a - b*sin(c + d*x)**4)
    F = -a**(sympy.S(3)/2)*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*b**(sympy.S(9)/4)*d*sqrt(sqrt(a) + sqrt(b))) - a**(sympy.S(3)/2)*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*b**(sympy.S(9)/4)*d*sqrt(sqrt(a) - sqrt(b))) + cos(c + d*x)**5/(5*b*d) - 2*cos(c + d*x)**3/(3*b*d) + (a + b)*cos(c + d*x)/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_195():
    f = sin(c + d*x)**7/(a - b*sin(c + d*x)**4)
    F = a*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*b**(sympy.S(7)/4)*d*sqrt(sqrt(a) + sqrt(b))) - a*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*b**(sympy.S(7)/4)*d*sqrt(sqrt(a) - sqrt(b))) - cos(c + d*x)**3/(3*b*d) + cos(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_196():
    f = sin(c + d*x)**5/(a - b*sin(c + d*x)**4)
    F = -sqrt(a)*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*b**(sympy.S(5)/4)*d*sqrt(sqrt(a) + sqrt(b))) - sqrt(a)*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*b**(sympy.S(5)/4)*d*sqrt(sqrt(a) - sqrt(b))) + cos(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_197():
    f = sin(c + d*x)**3/(a - b*sin(c + d*x)**4)
    F = atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*b**(sympy.S(3)/4)*d*sqrt(sqrt(a) + sqrt(b))) - atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*b**(sympy.S(3)/4)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_198():
    f = sin(c + d*x)/(a - b*sin(c + d*x)**4)
    F = -atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*sqrt(a)*b**(sympy.S(1)/4)*d*sqrt(sqrt(a) + sqrt(b))) - atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*sqrt(a)*b**(sympy.S(1)/4)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_199():
    f = csc(c + d*x)/(a - b*sin(c + d*x)**4)
    F = b**(sympy.S(1)/4)*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*a*d*sqrt(sqrt(a) + sqrt(b))) - b**(sympy.S(1)/4)*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*a*d*sqrt(sqrt(a) - sqrt(b))) - atanh(cos(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_200():
    f = csc(c + d*x)**3/(a - b*sin(c + d*x)**4)
    F = -atanh(cos(c + d*x))/(2*a*d) + 1/(4*a*d*(cos(c + d*x) + 1)) - 1/(4*a*d*(1 - cos(c + d*x))) - b**(sympy.S(3)/4)*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*a**(sympy.S(3)/2)*d*sqrt(sqrt(a) + sqrt(b))) - b**(sympy.S(3)/4)*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*a**(sympy.S(3)/2)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_201():
    f = csc(c + d*x)**5/(a - b*sin(c + d*x)**4)
    F = 3/(16*a*d*(cos(c + d*x) + 1)) + 1/(16*a*d*(cos(c + d*x) + 1)**2) - 3/(16*a*d*(1 - cos(c + d*x))) - 1/(16*a*d*(1 - cos(c + d*x))**2) + b**(sympy.S(5)/4)*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*a**2*d*sqrt(sqrt(a) + sqrt(b))) - b**(sympy.S(5)/4)*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*a**2*d*sqrt(sqrt(a) - sqrt(b))) - (3*a + 8*b)*atanh(cos(c + d*x))/(8*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_202():
    f = sin(c + d*x)**8/(a - b*sin(c + d*x)**4)
    F = a**(sympy.S(5)/4)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*b**2*d*sqrt(sqrt(a) + sqrt(b))) + a**(sympy.S(5)/4)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*b**2*d*sqrt(sqrt(a) - sqrt(b))) + 5*x/(8*b) - sin(c + d*x)*cos(c + d*x)**3/(4*b*d) + 5*sin(c + d*x)*cos(c + d*x)/(8*b*d) - x*(a + b)/b**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_203():
    f = sin(c + d*x)**6/(a - b*sin(c + d*x)**4)
    F = -a**(sympy.S(3)/4)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(3)/2)*d*sqrt(sqrt(a) + sqrt(b))) + a**(sympy.S(3)/4)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(3)/2)*d*sqrt(sqrt(a) - sqrt(b))) - x/(2*b) + sin(c + d*x)*cos(c + d*x)/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_204():
    f = sin(c + d*x)**4/(a - b*sin(c + d*x)**4)
    F = a**(sympy.S(1)/4)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*b*d*sqrt(sqrt(a) + sqrt(b))) + a**(sympy.S(1)/4)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*b*d*sqrt(sqrt(a) - sqrt(b))) - x/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_205():
    f = sin(c + d*x)**2/(a - b*sin(c + d*x)**4)
    F = -atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*sqrt(b)*d*sqrt(sqrt(a) + sqrt(b))) + atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*sqrt(b)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_206():
    f = 1/(a - b*sin(c + d*x)**4)
    F = atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*sqrt(sqrt(a) + sqrt(b))) + atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_207():
    f = csc(c + d*x)**2/(a - b*sin(c + d*x)**4)
    F = -cot(c + d*x)/(a*d) - sqrt(b)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)*d*sqrt(sqrt(a) + sqrt(b))) + sqrt(b)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_208():
    f = csc(c + d*x)**4/(a - b*sin(c + d*x)**4)
    F = -cot(c + d*x)**3/(3*a*d) - cot(c + d*x)/(a*d) + b*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4)*d*sqrt(sqrt(a) + sqrt(b))) + b*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_209():
    f = csc(c + d*x)**6/(a - b*sin(c + d*x)**4)
    F = -cot(c + d*x)**5/(5*a*d) - 2*cot(c + d*x)**3/(3*a*d) - (a + b)*cot(c + d*x)/(a**2*d) - b**(sympy.S(3)/2)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(9)/4)*d*sqrt(sqrt(a) + sqrt(b))) + b**(sympy.S(3)/2)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(9)/4)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_210():
    f = csc(c + d*x)**8/(a - b*sin(c + d*x)**4)
    F = -cot(c + d*x)**7/(7*a*d) - 3*cot(c + d*x)**5/(5*a*d) - (a + b)*cot(c + d*x)/(a**2*d) - (3*a + b)*cot(c + d*x)**3/(3*a**2*d) + b**2*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(11)/4)*d*sqrt(sqrt(a) + sqrt(b))) + b**2*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(11)/4)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_211():
    f = sin(c + d*x)**9/(a - b*sin(c + d*x)**4)**2
    F = sqrt(a)*(5*sqrt(a) + 6*sqrt(b))*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(8*b**(sympy.S(9)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) + sqrt(a)*(5*sqrt(a) - 6*sqrt(b))*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(8*b**(sympy.S(9)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2)) - a*(a - b*cos(c + d*x)**2 + b)*cos(c + d*x)/(b**2*d*(4*a - 4*b)*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)) - cos(c + d*x)/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_212():
    f = sin(c + d*x)**7/(a - b*sin(c + d*x)**4)**2
    F = -a*(2 - cos(c + d*x)**2)*cos(c + d*x)/(b*d*(4*a - 4*b)*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)) - (3*sqrt(a) + 4*sqrt(b))*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(8*b**(sympy.S(7)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) + (3*sqrt(a) - 4*sqrt(b))*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(8*b**(sympy.S(7)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_213():
    f = sin(c + d*x)**5/(a - b*sin(c + d*x)**4)**2
    F = -(a - b*cos(c + d*x)**2 + b)*cos(c + d*x)/(b*d*(4*a - 4*b)*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)) + (sqrt(a) - 2*sqrt(b))*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(8*sqrt(a)*b**(sympy.S(5)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2)) + (sqrt(a) + 2*sqrt(b))*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(8*sqrt(a)*b**(sympy.S(5)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_214():
    f = sin(c + d*x)**3/(a - b*sin(c + d*x)**4)**2
    F = -(2 - cos(c + d*x)**2)*cos(c + d*x)/(d*(4*a - 4*b)*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)) + atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(8*sqrt(a)*b**(sympy.S(3)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) - atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(8*sqrt(a)*b**(sympy.S(3)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_215():
    f = sin(c + d*x)/(a - b*sin(c + d*x)**4)**2
    F = -(a - b*cos(c + d*x)**2 + b)*cos(c + d*x)/(4*a*d*(a - b)*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)) - (3*sqrt(a) + 2*sqrt(b))*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(8*a**(sympy.S(3)/2)*b**(sympy.S(1)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) - (3*sqrt(a) - 2*sqrt(b))*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(8*a**(sympy.S(3)/2)*b**(sympy.S(1)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_216():
    f = csc(c + d*x)/(a - b*sin(c + d*x)**4)**2
    F = -b*(2 - cos(c + d*x)**2)*cos(c + d*x)/(4*a*d*(a - b)*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)) + b**(sympy.S(1)/4)*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*a**2*d*sqrt(sqrt(a) + sqrt(b))) - b**(sympy.S(1)/4)*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*a**2*d*sqrt(sqrt(a) - sqrt(b))) - atanh(cos(c + d*x))/(a**2*d) + b**(sympy.S(1)/4)*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(8*a**(sympy.S(3)/2)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) - b**(sympy.S(1)/4)*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(8*a**(sympy.S(3)/2)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_217():
    f = sin(c + d*x)**8/(a - b*sin(c + d*x)**4)**2
    F = -a**(sympy.S(1)/4)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*b**2*d*sqrt(sqrt(a) + sqrt(b))) - a**(sympy.S(1)/4)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*b**2*d*sqrt(sqrt(a) - sqrt(b))) - a**(sympy.S(1)/4)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(8*b**(sympy.S(3)/2)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) + a**(sympy.S(1)/4)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(8*b**(sympy.S(3)/2)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2)) + tan(c + d*x)**5/(4*b*d*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)) - tan(c + d*x)/(b*d*(4*a - 4*b)) + x/b**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_218():
    f = sin(c + d*x)**6/(a - b*sin(c + d*x)**4)**2
    F = tan(c + d*x)**3*sec(c + d*x)**2/(4*b*d*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)) - tan(c + d*x)/(b*d*(4*a - 4*b)) + (2*sqrt(a) + 3*sqrt(b))*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(1)/4)*b**(sympy.S(3)/2)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) - (2*sqrt(a) - 3*sqrt(b))*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(1)/4)*b**(sympy.S(3)/2)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_219():
    f = sin(c + d*x)**4/(a - b*sin(c + d*x)**4)**2
    F = tan(c + d*x)**5/(4*a*d*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)) - tan(c + d*x)/(4*a*d*(a - b)) - atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(3)/4)*sqrt(b)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) + atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(3)/4)*sqrt(b)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_220():
    f = sin(c + d*x)**2/(a - b*sin(c + d*x)**4)**2
    F = -(a + (a + b)*tan(c + d*x)**2)*tan(c + d*x)/(4*a*d*(a - b)*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)) - (2*sqrt(a) + sqrt(b))*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*sqrt(b)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) + (2*sqrt(a) - sqrt(b))*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*sqrt(b)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_221():
    f = (a - b*sin(c + d*x)**4)**(-2)
    F = -b*(2*tan(c + d*x)**2 + 1)*tan(c + d*x)/(4*a*d*(a - b)*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)) + (4*sqrt(a) + 3*sqrt(b))*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) + (4*sqrt(a) - 3*sqrt(b))*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_222():
    f = csc(c + d*x)**2/(a - b*sin(c + d*x)**4)**2
    F = -b*(a + (a + b)*tan(c + d*x)**2)*tan(c + d*x)/(4*a**2*d*(a - b)*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)) - cot(c + d*x)/(a**2*d) - sqrt(b)*(6*sqrt(a) + 5*sqrt(b))*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(9)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) + sqrt(b)*(6*sqrt(a) - 5*sqrt(b))*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(9)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_223():
    f = sin(c + d*x)**9/(a - b*sin(c + d*x)**4)**3
    F = -a*(a - b*cos(c + d*x)**2 + b)*cos(c + d*x)/(b**2*d*(8*a - 8*b)*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)**2) + (9*a**2 - 11*a*b - 10*b**2 - b*(4*a - 10*b)*cos(c + d*x)**2)*cos(c + d*x)/(32*b**2*d*(a - b)**2*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)) - (14*sqrt(a)*sqrt(b) + 5*a + 12*b)*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(64*sqrt(a)*b**(sympy.S(9)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) - (-14*sqrt(a)*sqrt(b) + 5*a + 12*b)*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(64*sqrt(a)*b**(sympy.S(9)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_224():
    f = sin(c + d*x)**7/(a - b*sin(c + d*x)**4)**3
    F = -a*(2 - cos(c + d*x)**2)*cos(c + d*x)/(b*d*(8*a - 8*b)*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)**2) + (5*a - 17*b - (3*a - 9*b)*cos(c + d*x)**2)*cos(c + d*x)/(32*b*d*(a - b)**2*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)) - (3*sqrt(a) + 6*sqrt(b))*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(64*sqrt(a)*b**(sympy.S(7)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) + (3*sqrt(a) - 6*sqrt(b))*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(64*sqrt(a)*b**(sympy.S(7)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_225():
    f = sin(c + d*x)**5/(a - b*sin(c + d*x)**4)**3
    F = -(a - b*cos(c + d*x)**2 + b)*cos(c + d*x)/(b*d*(8*a - 8*b)*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)**2) + (a**2 - 11*a*b - 2*b**2 + 2*b*(2*a + b)*cos(c + d*x)**2)*cos(c + d*x)/(32*a*b*d*(a - b)**2*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)) + (10*sqrt(a)*sqrt(b) + 3*a + 4*b)*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(64*a**(sympy.S(3)/2)*b**(sympy.S(5)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) + (-10*sqrt(a)*sqrt(b) + 3*a + 4*b)*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(64*a**(sympy.S(3)/2)*b**(sympy.S(5)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_226():
    f = sin(c + d*x)**3/(a - b*sin(c + d*x)**4)**3
    F = -(2 - cos(c + d*x)**2)*cos(c + d*x)/(d*(8*a - 8*b)*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)**2) - (11*a + b - (5*a + b)*cos(c + d*x)**2)*cos(c + d*x)/(32*a*d*(a - b)**2*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)) + (5*sqrt(a) + 2*sqrt(b))*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(64*a**(sympy.S(3)/2)*b**(sympy.S(3)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) - (5*sqrt(a) - 2*sqrt(b))*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(64*a**(sympy.S(3)/2)*b**(sympy.S(3)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_227():
    f = sin(c + d*x)/(a - b*sin(c + d*x)**4)**3
    F = -(a - b*cos(c + d*x)**2 + b)*cos(c + d*x)/(8*a*d*(a - b)*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)**2) - (-b*(12*a - 6*b)*cos(c + d*x)**2 + (a + 2*b)*(7*a - 3*b))*cos(c + d*x)/(32*a**2*d*(a - b)**2*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)) - (30*sqrt(a)*sqrt(b) + 21*a + 12*b)*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(64*a**(sympy.S(5)/2)*b**(sympy.S(1)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) - (-30*sqrt(a)*sqrt(b) + 21*a + 12*b)*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(64*a**(sympy.S(5)/2)*b**(sympy.S(1)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_228():
    f = csc(c + d*x)/(a - b*sin(c + d*x)**4)**3
    F = -b*(2 - cos(c + d*x)**2)*cos(c + d*x)/(8*a*d*(a - b)*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)**2) - b*(2 - cos(c + d*x)**2)*cos(c + d*x)/(4*a**2*d*(a - b)*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)) - b*(11*a + b - (5*a + b)*cos(c + d*x)**2)*cos(c + d*x)/(32*a**2*d*(a - b)**2*(a - b*cos(c + d*x)**4 + 2*b*cos(c + d*x)**2 - b)) + b**(sympy.S(1)/4)*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*a**3*d*sqrt(sqrt(a) + sqrt(b))) - b**(sympy.S(1)/4)*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*a**3*d*sqrt(sqrt(a) - sqrt(b))) - atanh(cos(c + d*x))/(a**3*d) + b**(sympy.S(1)/4)*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(8*a**(sympy.S(5)/2)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) + b**(sympy.S(1)/4)*(5*sqrt(a) + 2*sqrt(b))*atanh(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(64*a**(sympy.S(5)/2)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) - b**(sympy.S(1)/4)*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(8*a**(sympy.S(5)/2)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2)) - b**(sympy.S(1)/4)*(5*sqrt(a) - 2*sqrt(b))*atan(b**(sympy.S(1)/4)*cos(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(64*a**(sympy.S(5)/2)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_229():
    f = sin(c + d*x)**8/(a - b*sin(c + d*x)**4)**3
    F = tan(c + d*x)**9/(8*a*d*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)**2) - tan(c + d*x)**5*sec(c + d*x)**2/(32*a*b*d*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)) + tan(c + d*x)**3/(32*a*b*d*(a - b)) - (a + 5*b)*tan(c + d*x)/(32*a*b*d*(a - b)**2) + (2*sqrt(a) + 5*sqrt(b))*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(3)/4)*b**(sympy.S(3)/2)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) - (2*sqrt(a) - 5*sqrt(b))*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(3)/4)*b**(sympy.S(3)/2)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_230():
    f = sin(c + d*x)**6/(a - b*sin(c + d*x)**4)**3
    F = -(a*(a + 3*b) + (a**2 + 6*a*b + b**2)*tan(c + d*x)**2)*tan(c + d*x)/(8*d*(a - b)**3*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)**2) - (2*a*(a**2 - a*b - 8*b**2)/(a - b)**3 + (2*a**2 + 15*a*b + 3*b**2)*tan(c + d*x)**2/(a - b)**2)*tan(c + d*x)/(32*a*b*d*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)) + (10*sqrt(a)*sqrt(b) + 4*a + 3*b)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(5)/4)*b**(sympy.S(3)/2)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) - (-10*sqrt(a)*sqrt(b) + 4*a + 3*b)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(5)/4)*b**(sympy.S(3)/2)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_231():
    f = sin(c + d*x)**4/(a - b*sin(c + d*x)**4)**3
    F = -b*(3*a + b + (4*a + 4*b)*tan(c + d*x)**2)*tan(c + d*x)/(8*d*(a - b)**3*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)**2) - ((17*a + 3*b)*tan(c + d*x)**2/(a - b)**2 + (9*a**2 - 24*a*b - b**2)/(a - b)**3)*tan(c + d*x)/(32*a*d*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)) - (6*sqrt(a) + 3*sqrt(b))*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(7)/4)*sqrt(b)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) + (6*sqrt(a) - 3*sqrt(b))*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(7)/4)*sqrt(b)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_232():
    f = sin(c + d*x)**2/(a - b*sin(c + d*x)**4)**3
    F = -b*(a*(a + 3*b) + (a**2 + 6*a*b + b**2)*tan(c + d*x)**2)*tan(c + d*x)/(8*a*d*(a - b)**3*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)**2) - (2*a*(5*a**2 - 9*a*b - 4*b**2)/(a - b)**3 + (10*a**2 + 15*a*b - 5*b**2)*tan(c + d*x)**2/(a - b)**2)*tan(c + d*x)/(32*a**2*d*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)) - (14*sqrt(a)*sqrt(b) + 12*a + 5*b)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(9)/4)*sqrt(b)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) + (-14*sqrt(a)*sqrt(b) + 12*a + 5*b)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(9)/4)*sqrt(b)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_233():
    f = (a - b*sin(c + d*x)**4)**(-3)
    F = -b**2*(3*a + b + (4*a + 4*b)*tan(c + d*x)**2)*tan(c + d*x)/(8*a*d*(a - b)**3*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)**2) - b*((33*a - 13*b)*tan(c + d*x)**2/(a - b)**2 + (17*a**2 - 40*a*b + 7*b**2)/(a - b)**3)*tan(c + d*x)/(32*a**2*d*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)) + (50*sqrt(a)*sqrt(b) + 32*a + 21*b)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) + (-50*sqrt(a)*sqrt(b) + 32*a + 21*b)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_234():
    f = csc(c + d*x)**2/(a - b*sin(c + d*x)**4)**3
    F = -b**2*(a*(a + 3*b) + (a**2 + 6*a*b + b**2)*tan(c + d*x)**2)*tan(c + d*x)/(8*a**2*d*(a - b)**3*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)**2) - b*(2*a**2*(9*a - 17*b)/(a - b)**3 + (18*a**2 + 15*a*b - 13*b**2)*tan(c + d*x)**2/(a - b)**2)*tan(c + d*x)/(32*a**3*d*(2*a*tan(c + d*x)**2 + a + (a - b)*tan(c + d*x)**4)) - cot(c + d*x)/(a**3*d) - 3*sqrt(b)*(34*sqrt(a)*sqrt(b) + 20*a + 15*b)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(13)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) + 3*sqrt(b)*(-34*sqrt(a)*sqrt(b) + 20*a + 15*b)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(13)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_235():
    f = 1/(a + b*sin(x)**4)
    F = sqrt(2)*(sqrt(a) - sqrt(a + b))*log(-sqrt(2)*a**(sympy.S(1)/4)*sqrt(-sqrt(a)*sqrt(a + b) + a + b)*tan(x) + sqrt(a)*(a + b)**(sympy.S(1)/4) + (a + b)**(sympy.S(3)/4)*tan(x)**2)/(8*a**(sympy.S(3)/4)*(a + b)**(sympy.S(1)/4)*sqrt(-sqrt(a)*sqrt(a + b) + a + b)) - sqrt(2)*(sqrt(a) - sqrt(a + b))*log(sqrt(2)*a**(sympy.S(1)/4)*sqrt(-sqrt(a)*sqrt(a + b) + a + b)*tan(x) + sqrt(a)*(a + b)**(sympy.S(1)/4) + (a + b)**(sympy.S(3)/4)*tan(x)**2)/(8*a**(sympy.S(3)/4)*(a + b)**(sympy.S(1)/4)*sqrt(-sqrt(a)*sqrt(a + b) + a + b)) - sqrt(2)*(sqrt(a) + sqrt(a + b))*atan((a**(sympy.S(1)/4)*sqrt(-sqrt(a)*sqrt(a + b) + a + b) - sqrt(2)*(a + b)**(sympy.S(3)/4)*tan(x))/(a**(sympy.S(1)/4)*sqrt(sqrt(a)*sqrt(a + b) + a + b)))/(4*a**(sympy.S(3)/4)*(a + b)**(sympy.S(1)/4)*sqrt(sqrt(a)*sqrt(a + b) + a + b)) + sqrt(2)*(sqrt(a) + sqrt(a + b))*atan((a**(sympy.S(1)/4)*sqrt(-sqrt(a)*sqrt(a + b) + a + b) + sqrt(2)*(a + b)**(sympy.S(3)/4)*tan(x))/(a**(sympy.S(1)/4)*sqrt(sqrt(a)*sqrt(a + b) + a + b)))/(4*a**(sympy.S(3)/4)*(a + b)**(sympy.S(1)/4)*sqrt(sqrt(a)*sqrt(a + b) + a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_236():
    f = 1/(sin(x)**4 + 1)
    F = x/(2*sqrt(-1 + sqrt(2))) + sqrt(-1 + sqrt(2))*log(sqrt(2)*tan(x)**2 + sqrt(-2 + 2*sqrt(2))*tan(x) + 1)/8 - sqrt(-1 + sqrt(2))*log(2*tan(x)**2 - 2*sqrt(-1 + sqrt(2))*tan(x) + sqrt(2))/8 + atan((-(-2 + sqrt(2))*sin(x)*cos(x) - 2*sqrt(-1 + sqrt(2))*cos(x)**2 + sqrt(-1 + sqrt(2)))/(-2*sqrt(-1 + sqrt(2))*sin(x)*cos(x) + (-2 + sqrt(2))*cos(x)**2 + sqrt(1 + sqrt(2)) + 2))/(4*sqrt(-1 + sqrt(2))) - atan(((-2 + sqrt(2))*sin(x)*cos(x) - 2*sqrt(-1 + sqrt(2))*cos(x)**2 + sqrt(-1 + sqrt(2)))/(2*sqrt(-1 + sqrt(2))*sin(x)*cos(x) + (-2 + sqrt(2))*cos(x)**2 + sqrt(1 + sqrt(2)) + 2))/(4*sqrt(-1 + sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_237():
    f = sqrt(a + b*sin(c + d*x)**4)*sin(c + d*x)
    F = -2*b**(sympy.S(1)/4)*sqrt((a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b)/((a + b)*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)**2))*(a + b)**(sympy.S(3)/4)*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)*elliptic_e(2*atan(b**(sympy.S(1)/4)*cos(c + d*x)/(a + b)**(sympy.S(1)/4)), sqrt(b)/(2*sqrt(a + b)) + sympy.S.Half)/(3*d*sqrt(a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b)) + 2*sqrt(b)*sqrt(a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b)*cos(c + d*x)/(3*d*sqrt(a + b)*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)) - sqrt(a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b)*cos(c + d*x)/(3*d) + sqrt((a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b)/((a + b)*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)**2))*(a + b)**(sympy.S(3)/4)*(sqrt(b) - sqrt(a + b))*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)*elliptic_f(2*atan(b**(sympy.S(1)/4)*cos(c + d*x)/(a + b)**(sympy.S(1)/4)), sqrt(b)/(2*sqrt(a + b)) + sympy.S.Half)/(3*b**(sympy.S(1)/4)*d*sqrt(a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_238():
    f = sqrt(a + b*sin(c + d*x)**4)*csc(c + d*x)
    F = ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.atan(((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4))))))**(Integer(-1))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt(Symbol('b')) * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * (Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(4))**(Integer(-1))) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))) * (((Symbol('a') + Symbol('b')) * ((Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))))**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(2))**(Integer(-1)) * (Integer(1) + (sympy.sqrt(Symbol('b')) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1))))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('b')) + (Integer(-1) * sympy.sqrt((Symbol('a') + Symbol('b'))))))**(Integer(2)) * (Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))) * (((Symbol('a') + Symbol('b')) * ((Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))))**(Integer(2))))**(Integer(-1)))) * sympy.Function('EllipticPi')((((sympy.sqrt(Symbol('b')) + sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(2))**(Integer(-1)) * (Integer(1) + (sympy.sqrt(Symbol('b')) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1))))))) * ((Integer(4) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_239():
    f = sin(c + d*x)**5/sqrt(a + b*sin(c + d*x)**4)
    F = -sqrt(a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b)*cos(c + d*x)/(3*b*d) + 2*sqrt(a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b)*cos(c + d*x)/(3*sqrt(b)*d*sqrt(a + b)*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)) - 2*sqrt((a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b)/((a + b)*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)**2))*(a + b)**(sympy.S(3)/4)*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)*elliptic_e(2*atan(b**(sympy.S(1)/4)*cos(c + d*x)/(a + b)**(sympy.S(1)/4)), sqrt(b)/(2*sqrt(a + b)) + sympy.S.Half)/(3*b**(sympy.S(3)/4)*d*sqrt(a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b)) + sqrt((a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b)/((a + b)*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)**2))*(a + b)**(sympy.S(1)/4)*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)*(a + 2*sqrt(b)*sqrt(a + b) - 2*b)*elliptic_f(2*atan(b**(sympy.S(1)/4)*cos(c + d*x)/(a + b)**(sympy.S(1)/4)), sqrt(b)/(2*sqrt(a + b)) + sympy.S.Half)/(6*b**(sympy.S(5)/4)*d*sqrt(a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_240():
    f = sin(c + d*x)**3/sqrt(a + b*sin(c + d*x)**4)
    F = sqrt(a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b)*cos(c + d*x)/(sqrt(b)*d*sqrt(a + b)*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)) - sqrt((a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b)/((a + b)*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)**2))*(a + b)**(sympy.S(3)/4)*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)*elliptic_e(2*atan(b**(sympy.S(1)/4)*cos(c + d*x)/(a + b)**(sympy.S(1)/4)), sqrt(b)/(2*sqrt(a + b)) + sympy.S.Half)/(b**(sympy.S(3)/4)*d*sqrt(a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b)) - sqrt((a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b)/((a + b)*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)**2))*(a + b)**(sympy.S(1)/4)*(sqrt(b) - sqrt(a + b))*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)*elliptic_f(2*atan(b**(sympy.S(1)/4)*cos(c + d*x)/(a + b)**(sympy.S(1)/4)), sqrt(b)/(2*sqrt(a + b)) + sympy.S.Half)/(2*b**(sympy.S(3)/4)*d*sqrt(a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_241():
    f = sin(c + d*x)/sqrt(a + b*sin(c + d*x)**4)
    F = -sqrt((a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b)/((a + b)*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)**2))*(a + b)**(sympy.S(1)/4)*(sqrt(b)*cos(c + d*x)**2/sqrt(a + b) + 1)*elliptic_f(2*atan(b**(sympy.S(1)/4)*cos(c + d*x)/(a + b)**(sympy.S(1)/4)), sqrt(b)/(2*sqrt(a + b)) + sympy.S.Half)/(2*b**(sympy.S(1)/4)*d*sqrt(a + b*cos(c + d*x)**4 - 2*b*cos(c + d*x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_242():
    f = csc(c + d*x)/sqrt(a + b*sin(c + d*x)**4)
    F = (Integer(-1) * (sympy.atan(((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4))))))**(Integer(-1)))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('d')))**(Integer(-1)))) + (((Symbol('b'))**((Integer(4))**(Integer(-1))) * ((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('b')) + (Integer(-1) * sympy.sqrt((Symbol('a') + Symbol('b'))))) * (Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))) * (((Symbol('a') + Symbol('b')) * ((Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))))**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(2))**(Integer(-1)) * (Integer(1) + (sympy.sqrt(Symbol('b')) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1))))))) * ((Integer(2) * Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('b')) + (Integer(-1) * sympy.sqrt((Symbol('a') + Symbol('b'))))))**(Integer(2)) * (Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))) * (((Symbol('a') + Symbol('b')) * ((Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))))**(Integer(2))))**(Integer(-1)))) * sympy.Function('EllipticPi')((((sympy.sqrt(Symbol('b')) + sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(2))**(Integer(-1)) * (Integer(1) + (sympy.sqrt(Symbol('b')) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1))))))) * ((Integer(4) * Symbol('a') * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_243():
    f = csc(c + d*x)**3/sqrt(a + b*sin(c + d*x)**4)
    F = (Integer(-1) * (sympy.atan(((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4))))))**(Integer(-1)))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * (Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.csc((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * Symbol('d')))**(Integer(-1)))) + (((Symbol('b'))**((Integer(4))**(Integer(-1))) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))) * (((Symbol('a') + Symbol('b')) * ((Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))))**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(2))**(Integer(-1)) * (Integer(1) + (sympy.sqrt(Symbol('b')) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1))))))) * ((Integer(2) * Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('a') + Symbol('b') + (Integer(-1) * (sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('a') + Symbol('b')))))) * (Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))) * (((Symbol('a') + Symbol('b')) * ((Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))))**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(2))**(Integer(-1)) * (Integer(1) + (sympy.sqrt(Symbol('b')) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1))))))) * ((Integer(2) * Symbol('a') * ((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('b')) + (Integer(-1) * sympy.sqrt((Symbol('a') + Symbol('b'))))))**(Integer(2)) * (Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))) * (((Symbol('a') + Symbol('b')) * ((Integer(1) + ((sympy.sqrt(Symbol('b')) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))))**(Integer(2))))**(Integer(-1)))) * sympy.Function('EllipticPi')((((sympy.sqrt(Symbol('b')) + sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(2))**(Integer(-1)) * (Integer(1) + (sympy.sqrt(Symbol('b')) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1))))))) * ((Integer(8) * Symbol('a') * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt((Symbol('a') + Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) + (Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_244():
    f = sin(c + d*x)**2/sqrt(a + b*sin(c + d*x)**4)
    F = (Integer(-1) * ((sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt((Symbol('a') + (Integer(2) * Symbol('a') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) + ((Symbol('a') + Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4))))))**(Integer(-1)))) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(2) * Symbol('a') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) + ((Symbol('a') + Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + sympy.sqrt((Symbol('a') + Symbol('b')))) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.elliptic_f((Integer(2) * sympy.atan(((((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(2))**(Integer(-1)) * (Integer(1) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1))))))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt((Symbol('a') + Symbol('b'))) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Integer(2) * Symbol('a') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) + ((Symbol('a') + Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt((Symbol('a') + Symbol('b'))) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * ((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1)))) + ((((sympy.sqrt(Symbol('a')) + sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(2)) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.Function('EllipticPi')((Integer(-1) * (((sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt((Symbol('a') + Symbol('b'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1)))), (Integer(2) * sympy.atan(((((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(2))**(Integer(-1)) * (Integer(1) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1))))))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt((Symbol('a') + Symbol('b'))) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Integer(2) * Symbol('a') * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) + ((Symbol('a') + Symbol('b')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt((Symbol('a') + Symbol('b'))) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_245():
    f = 1/sqrt(a + b*sin(c + d*x)**4)
    F = sqrt((2*a*tan(c + d*x)**2 + a + (a + b)*tan(c + d*x)**4)/(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)**2)*(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)*cos(c + d*x)**2*elliptic_f(2*atan((a + b)**(sympy.S(1)/4)*tan(c + d*x)/a**(sympy.S(1)/4)), -sqrt(a)/(2*sqrt(a + b)) + sympy.S.Half)/(2*a**(sympy.S(1)/4)*d*(a + b)**(sympy.S(1)/4)*sqrt(a + b*sin(c + d*x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_246():
    f = csc(c + d*x)**2/sqrt(a + b*sin(c + d*x)**4)
    F = -(2*a*tan(c + d*x)**2 + a + (a + b)*tan(c + d*x)**4)*cos(c + d*x)**2*cot(c + d*x)/(a*d*sqrt(a + b*sin(c + d*x)**4)) + sqrt(a + b)*(2*a*tan(c + d*x)**2 + a + (a + b)*tan(c + d*x)**4)*sin(c + d*x)*cos(c + d*x)/(a*d*(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)*sqrt(a + b*sin(c + d*x)**4)) - sqrt((2*a*tan(c + d*x)**2 + a + (a + b)*tan(c + d*x)**4)/(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)**2)*(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)*(a + b)**(sympy.S(1)/4)*cos(c + d*x)**2*elliptic_e(2*atan((a + b)**(sympy.S(1)/4)*tan(c + d*x)/a**(sympy.S(1)/4)), -sqrt(a)/(2*sqrt(a + b)) + sympy.S.Half)/(a**(sympy.S(3)/4)*d*sqrt(a + b*sin(c + d*x)**4)) + sqrt((2*a*tan(c + d*x)**2 + a + (a + b)*tan(c + d*x)**4)/(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)**2)*(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)*(sqrt(a)*sqrt(a + b) + a + b)*cos(c + d*x)**2*elliptic_f(2*atan((a + b)**(sympy.S(1)/4)*tan(c + d*x)/a**(sympy.S(1)/4)), -sqrt(a)/(2*sqrt(a + b)) + sympy.S.Half)/(2*a**(sympy.S(3)/4)*d*(a + b)**(sympy.S(3)/4)*sqrt(a + b*sin(c + d*x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_247():
    f = 1/(a + b*sin(x)**5)
    F = 2*atan((a**(sympy.S(1)/5)*tan(x/2) + (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5))/sqrt(a**(sympy.S(2)/5) - (-1)**(sympy.S(4)/5)*b**(sympy.S(2)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(2)/5) - (-1)**(sympy.S(4)/5)*b**(sympy.S(2)/5))) + 2*atan((a**(sympy.S(1)/5)*tan(x/2) + (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5))/sqrt(a**(sympy.S(2)/5) + (-1)**(sympy.S(3)/5)*b**(sympy.S(2)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(2)/5) + (-1)**(sympy.S(3)/5)*b**(sympy.S(2)/5))) - 2*atan((-1)**(sympy.S(1)/5)*((-1)**(sympy.S(4)/5)*a**(sympy.S(1)/5)*tan(x/2) + b**(sympy.S(1)/5))/sqrt(a**(sympy.S(2)/5) - (-1)**(sympy.S(2)/5)*b**(sympy.S(2)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(2)/5) - (-1)**(sympy.S(2)/5)*b**(sympy.S(2)/5))) - 2*I*atanh((-1)**(sympy.S(1)/10)*((-1)**(sympy.S(2)/5)*a**(sympy.S(1)/5)*tan(x/2) + b**(sympy.S(1)/5))/sqrt(a**(sympy.S(2)/5) + (-1)**(sympy.S(1)/5)*b**(sympy.S(2)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(2)/5) + (-1)**(sympy.S(1)/5)*b**(sympy.S(2)/5))) + 2*atan((a**(sympy.S(1)/5)*tan(x/2) + b**(sympy.S(1)/5))/sqrt(a**(sympy.S(2)/5) - b**(sympy.S(2)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(2)/5) - b**(sympy.S(2)/5)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_248():
    f = 1/(a + b*sin(x)**6)
    F = atan(sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))*tan(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))) + atan(sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))*tan(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))) + atan(sqrt(a**(sympy.S(1)/3) + b**(sympy.S(1)/3))*tan(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_249():
    f = 1/(a + b*sin(x)**8)
    F = -atan(sqrt(I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))*tan(x)/(-a)**(sympy.S(1)/8))/(4*(-a)**(sympy.S(7)/8)*sqrt(I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))) - atan(sqrt(-I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))*tan(x)/(-a)**(sympy.S(1)/8))/(4*(-a)**(sympy.S(7)/8)*sqrt(-I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))) - atan(sqrt(b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))*tan(x)/(-a)**(sympy.S(1)/8))/(4*(-a)**(sympy.S(7)/8)*sqrt(b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))) - atan(sqrt(-b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))*tan(x)/(-a)**(sympy.S(1)/8))/(4*(-a)**(sympy.S(7)/8)*sqrt(-b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_250():
    f = 1/(a - b*sin(x)**5)
    F = -2*atan((-a**(sympy.S(1)/5)*tan(x/2) + (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5))/sqrt(a**(sympy.S(2)/5) - (-1)**(sympy.S(4)/5)*b**(sympy.S(2)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(2)/5) - (-1)**(sympy.S(4)/5)*b**(sympy.S(2)/5))) - 2*atan((-a**(sympy.S(1)/5)*tan(x/2) + (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5))/sqrt(a**(sympy.S(2)/5) + (-1)**(sympy.S(3)/5)*b**(sympy.S(2)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(2)/5) + (-1)**(sympy.S(3)/5)*b**(sympy.S(2)/5))) + 2*atan((a**(sympy.S(1)/5)*tan(x/2) + (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5))/sqrt(a**(sympy.S(2)/5) - (-1)**(sympy.S(2)/5)*b**(sympy.S(2)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(2)/5) - (-1)**(sympy.S(2)/5)*b**(sympy.S(2)/5))) + 2*atan((a**(sympy.S(1)/5)*tan(x/2) + (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5))/sqrt(a**(sympy.S(2)/5) + (-1)**(sympy.S(1)/5)*b**(sympy.S(2)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(2)/5) + (-1)**(sympy.S(1)/5)*b**(sympy.S(2)/5))) - 2*atan((-a**(sympy.S(1)/5)*tan(x/2) + b**(sympy.S(1)/5))/sqrt(a**(sympy.S(2)/5) - b**(sympy.S(2)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(2)/5) - b**(sympy.S(2)/5)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_251():
    f = 1/(a - b*sin(x)**6)
    F = atan(sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))*tan(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))) + atan(sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))*tan(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))) + atan(sqrt(a**(sympy.S(1)/3) - b**(sympy.S(1)/3))*tan(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_252():
    f = 1/(a - b*sin(x)**8)
    F = atan(sqrt(a**(sympy.S(1)/4) + I*b**(sympy.S(1)/4))*tan(x)/a**(sympy.S(1)/8))/(4*a**(sympy.S(7)/8)*sqrt(a**(sympy.S(1)/4) + I*b**(sympy.S(1)/4))) + atan(sqrt(a**(sympy.S(1)/4) - I*b**(sympy.S(1)/4))*tan(x)/a**(sympy.S(1)/8))/(4*a**(sympy.S(7)/8)*sqrt(a**(sympy.S(1)/4) - I*b**(sympy.S(1)/4))) + atan(sqrt(a**(sympy.S(1)/4) + b**(sympy.S(1)/4))*tan(x)/a**(sympy.S(1)/8))/(4*a**(sympy.S(7)/8)*sqrt(a**(sympy.S(1)/4) + b**(sympy.S(1)/4))) + atan(sqrt(a**(sympy.S(1)/4) - b**(sympy.S(1)/4))*tan(x)/a**(sympy.S(1)/8))/(4*a**(sympy.S(7)/8)*sqrt(a**(sympy.S(1)/4) - b**(sympy.S(1)/4)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_253():
    f = 1/(sin(x)**5 + 1)
    F = 2*atan((tan(x/2) + (-1)**(sympy.S(4)/5))/sqrt(1 + (-1)**(sympy.S(3)/5)))/(5*sqrt(1 + (-1)**(sympy.S(3)/5))) + 2*atan((tan(x/2) + (-1)**(sympy.S(2)/5))/sqrt(1 - (-1)**(sympy.S(4)/5)))/(5*sqrt(1 - (-1)**(sympy.S(4)/5))) - 2*atan((-1)**(sympy.S(1)/5)*((-1)**(sympy.S(4)/5)*tan(x/2) + 1)/sqrt(1 - (-1)**(sympy.S(2)/5)))/(5*sqrt(1 - (-1)**(sympy.S(2)/5))) - 2*I*atanh((-1)**(sympy.S(1)/10)*((-1)**(sympy.S(2)/5)*tan(x/2) + 1)/sqrt(1 + (-1)**(sympy.S(1)/5)))/(5*sqrt(1 + (-1)**(sympy.S(1)/5))) - cos(x)/(5*sin(x) + 5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_254():
    f = 1/(sin(x)**6 + 1)
    F = sqrt(2)*x/6 + atan(sqrt(1 - (-1)**(sympy.S(1)/3))*tan(x))/(3*sqrt(1 - (-1)**(sympy.S(1)/3))) + atan(sqrt(1 + (-1)**(sympy.S(2)/3))*tan(x))/(3*sqrt(1 + (-1)**(sympy.S(2)/3))) + sqrt(2)*atan(sin(x)*cos(x)/(sin(x)**2 + 1 + sqrt(2)))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_255():
    f = 1/(1 - sin(x)**5)
    F = 2*atan((tan(x/2) + (-1)**(sympy.S(3)/5))/sqrt(1 + (-1)**(sympy.S(1)/5)))/(5*sqrt(1 + (-1)**(sympy.S(1)/5))) + 2*atan((tan(x/2) + (-1)**(sympy.S(1)/5))/sqrt(1 - (-1)**(sympy.S(2)/5)))/(5*sqrt(1 - (-1)**(sympy.S(2)/5))) - 2*atan((-tan(x/2) + (-1)**(sympy.S(4)/5))/sqrt(1 + (-1)**(sympy.S(3)/5)))/(5*sqrt(1 + (-1)**(sympy.S(3)/5))) - 2*atan((-tan(x/2) + (-1)**(sympy.S(2)/5))/sqrt(1 - (-1)**(sympy.S(4)/5)))/(5*sqrt(1 - (-1)**(sympy.S(4)/5))) + cos(x)/(5 - 5*sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_256():
    f = 1/(1 - sin(x)**6)
    F = tan(x)/3 + atan(sqrt(1 + (-1)**(sympy.S(1)/3))*tan(x))/(3*sqrt(1 + (-1)**(sympy.S(1)/3))) + atan(sqrt(1 - (-1)**(sympy.S(2)/3))*tan(x))/(3*sqrt(1 - (-1)**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_257():
    f = 1/(1 - sin(x)**8)
    F = sqrt(2)*x/8 + tan(x)/4 + atan(sqrt(1 - I)*tan(x))/(4*sqrt(1 - I)) + atan(sqrt(1 + I)*tan(x))/(4*sqrt(1 + I)) + sqrt(2)*atan(sin(x)*cos(x)/(sin(x)**2 + 1 + sqrt(2)))/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_258():
    f = cos(x)**9/(-a*sin(x)**2 + a)
    F = -sin(x)**7/(7*a) + 3*sin(x)**5/(5*a) - sin(x)**3/a + sin(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_259():
    f = cos(x)**7/(-a*sin(x)**2 + a)
    F = sin(x)**5/(5*a) - 2*sin(x)**3/(3*a) + sin(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_260():
    f = cos(x)**5/(-a*sin(x)**2 + a)
    F = -sin(x)**3/(3*a) + sin(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_261():
    f = cos(x)**3/(-a*sin(x)**2 + a)
    F = sin(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_262():
    f = cos(x)/(-a*sin(x)**2 + a)
    F = atanh(sin(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_263():
    f = sec(x)**3/(-a*sin(x)**2 + a)
    F = tan(x)*sec(x)**3/(4*a) + 3*tan(x)*sec(x)/(8*a) + 3*atanh(sin(x))/(8*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_264():
    f = cos(x)**6/(-a*sin(x)**2 + a)
    F = 3*x/(8*a) + sin(x)*cos(x)**3/(4*a) + 3*sin(x)*cos(x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_265():
    f = cos(x)**4/(-a*sin(x)**2 + a)
    F = x/(2*a) + sin(x)*cos(x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_266():
    f = cos(x)**2/(-a*sin(x)**2 + a)
    F = x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_267():
    f = sec(x)/(-a*sin(x)**2 + a)
    F = tan(x)*sec(x)/(2*a) + atanh(sin(x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_268():
    f = sec(x)**2/(-a*sin(x)**2 + a)
    F = tan(x)**3/(3*a) + tan(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_269():
    f = sec(x)**4/(-a*sin(x)**2 + a)
    F = tan(x)**5/(5*a) + 2*tan(x)**3/(3*a) + tan(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_270():
    f = cos(x)**9/(-a*sin(x)**2 + a)**2
    F = sin(x)**5/(5*a**2) - 2*sin(x)**3/(3*a**2) + sin(x)/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_271():
    f = cos(x)**7/(-a*sin(x)**2 + a)**2
    F = -sin(x)**3/(3*a**2) + sin(x)/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_272():
    f = cos(x)**5/(-a*sin(x)**2 + a)**2
    F = sin(x)/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_273():
    f = cos(x)**3/(-a*sin(x)**2 + a)**2
    F = atanh(sin(x))/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_274():
    f = cos(x)/(-a*sin(x)**2 + a)**2
    F = tan(x)*sec(x)/(2*a**2) + atanh(sin(x))/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_275():
    f = sec(x)/(-a*sin(x)**2 + a)**2
    F = tan(x)*sec(x)**3/(4*a**2) + 3*tan(x)*sec(x)/(8*a**2) + 3*atanh(sin(x))/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_276():
    f = cos(x)**8/(-a*sin(x)**2 + a)**2
    F = 3*x/(8*a**2) + sin(x)*cos(x)**3/(4*a**2) + 3*sin(x)*cos(x)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_277():
    f = cos(x)**6/(-a*sin(x)**2 + a)**2
    F = x/(2*a**2) + sin(x)*cos(x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_278():
    f = cos(x)**4/(-a*sin(x)**2 + a)**2
    F = x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_279():
    f = cos(x)**2/(-a*sin(x)**2 + a)**2
    F = tan(x)/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_280():
    f = sec(x)**2/(-a*sin(x)**2 + a)**2
    F = tan(x)**5/(5*a**2) + 2*tan(x)**3/(3*a**2) + tan(x)/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_281():
    f = sec(x)**4/(-a*sin(x)**2 + a)**2
    F = tan(x)**7/(7*a**2) + 3*tan(x)**5/(5*a**2) + tan(x)**3/a**2 + tan(x)/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_282():
    f = (a + b*sin(e + f*x)**2)*cos(e + f*x)**6
    F = -b*sin(e + f*x)*cos(e + f*x)**7/(8*f) + x*(5*a/16 + 5*b/128) + (8*a + b)*sin(e + f*x)*cos(e + f*x)**5/(48*f) + (40*a + 5*b)*sin(e + f*x)*cos(e + f*x)**3/(192*f) + (40*a + 5*b)*sin(e + f*x)*cos(e + f*x)/(128*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_283():
    f = (a + b*sin(e + f*x)**2)*cos(e + f*x)**4
    F = -b*sin(e + f*x)*cos(e + f*x)**5/(6*f) + x*(3*a/8 + b/16) + (6*a + b)*sin(e + f*x)*cos(e + f*x)**3/(24*f) + (6*a + b)*sin(e + f*x)*cos(e + f*x)/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_284():
    f = (a + b*sin(e + f*x)**2)*cos(e + f*x)**2
    F = -b*sin(e + f*x)*cos(e + f*x)**3/(4*f) + x*(a/2 + b/8) + (4*a + b)*sin(e + f*x)*cos(e + f*x)/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_285():
    f = a + b*sin(e + f*x)**2
    F = a*x + b*x/2 - b*sin(e + f*x)*cos(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_286():
    f = (a + b*sin(e + f*x)**2)*sec(e + f*x)**2
    F = -b*x + (a + b)*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_287():
    f = (a + b*sin(e + f*x)**2)*sec(e + f*x)**4
    F = a*tan(e + f*x)/f + (a + b)*tan(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_288():
    f = (a + b*sin(e + f*x)**2)*sec(e + f*x)**6
    F = a*tan(e + f*x)/f + (a + b)*tan(e + f*x)**5/(5*f) + (2*a + b)*tan(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_289():
    f = (a + b*sin(e + f*x)**2)*sec(e + f*x)**8
    F = a*tan(e + f*x)/f + (a + b)*tan(e + f*x)**7/(7*f) + (3*a + b)*tan(e + f*x)**3/(3*f) + (3*a + 2*b)*tan(e + f*x)**5/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_290():
    f = (a + b*sin(e + f*x)**2)**2*cos(e + f*x)**4
    F = -b*(a + (a + b)*tan(e + f*x)**2)*sin(e + f*x)*cos(e + f*x)**7/(8*f) - b*(10*a + 3*b)*sin(e + f*x)*cos(e + f*x)**5/(48*f) + x*(3*a**2/8 + a*b/8 + 3*b**2/128) + (48*a**2 + 16*a*b + 3*b**2)*sin(e + f*x)*cos(e + f*x)**3/(192*f) + (48*a**2 + 16*a*b + 3*b**2)*sin(e + f*x)*cos(e + f*x)/(128*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_291():
    f = (a + b*sin(e + f*x)**2)**2*cos(e + f*x)**2
    F = -b*(a + (a + b)*tan(e + f*x)**2)*sin(e + f*x)*cos(e + f*x)**5/(6*f) - b*(8*a + 3*b)*sin(e + f*x)*cos(e + f*x)**3/(24*f) + x*(a**2/2 + a*b/4 + b**2/16) + (8*a**2 + 4*a*b + b**2)*sin(e + f*x)*cos(e + f*x)/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_292():
    f = (a + b*sin(e + f*x)**2)**2
    F = -b**2*sin(e + f*x)**3*cos(e + f*x)/(4*f) - b*(8*a + 3*b)*sin(e + f*x)*cos(e + f*x)/(8*f) + x*(a**2 + a*b + 3*b**2/8)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_293():
    f = (a + b*sin(e + f*x)**2)**2*sec(e + f*x)**2
    F = b**2*sin(e + f*x)*cos(e + f*x)/(2*f) - b*x*(4*a + 3*b)/2 + (a + b)**2*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_294():
    f = (a + b*sin(e + f*x)**2)**2*sec(e + f*x)**4
    F = b**2*x + (a + b)**2*tan(e + f*x)**3/(3*f) + (a**2 - b**2)*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_295():
    f = (a + b*sin(e + f*x)**2)**2*sec(e + f*x)**6
    F = a**2*tan(e + f*x)/f + 2*a*(a + b)*tan(e + f*x)**3/(3*f) + (a + b)**2*tan(e + f*x)**5/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_296():
    f = (a + b*sin(e + f*x)**2)**2*sec(e + f*x)**8
    F = a**2*tan(e + f*x)/f + a*(3*a + 2*b)*tan(e + f*x)**3/(3*f) + (a + b)**2*tan(e + f*x)**7/(7*f) + (a + b)*(3*a + b)*tan(e + f*x)**5/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_297():
    f = (a + b*sin(e + f*x)**2)**2*sec(e + f*x)**10
    F = a**2*tan(e + f*x)/f + 2*a*(2*a + b)*tan(e + f*x)**3/(3*f) + (a + b)**2*tan(e + f*x)**9/(9*f) + (2*a + b)*(2*a + 2*b)*tan(e + f*x)**7/(7*f) + (6*a**2 + 6*a*b + b**2)*tan(e + f*x)**5/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_298():
    f = cos(x)**7/(a + b*sin(x)**2)
    F = -sin(x)**5/(5*b) + (a + 3*b)*sin(x)**3/(3*b**2) - (a**2 + 3*a*b + 3*b**2)*sin(x)/b**3 + (a + b)**3*atan(sqrt(b)*sin(x)/sqrt(a))/(sqrt(a)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_299():
    f = cos(x)**6/(a + b*sin(x)**2)
    F = -sin(x)*cos(x)**3/(4*b) - (4*a + 7*b)*sin(x)*cos(x)/(8*b**2) - x*(8*a**2 + 20*a*b + 15*b**2)/(8*b**3) + (a + b)**(sympy.S(5)/2)*atan(sqrt(a + b)*tan(x)/sqrt(a))/(sqrt(a)*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_300():
    f = cos(x)**5/(a + b*sin(x)**2)
    F = sin(x)**3/(3*b) - (a + 2*b)*sin(x)/b**2 + (a + b)**2*atan(sqrt(b)*sin(x)/sqrt(a))/(sqrt(a)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_301():
    f = cos(x)**4/(a + b*sin(x)**2)
    F = -sin(x)*cos(x)/(2*b) - x*(2*a + 3*b)/(2*b**2) + (a + b)**(sympy.S(3)/2)*atan(sqrt(a + b)*tan(x)/sqrt(a))/(sqrt(a)*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_302():
    f = cos(x)**3/(a + b*sin(x)**2)
    F = -sin(x)/b + (a + b)*atan(sqrt(b)*sin(x)/sqrt(a))/(sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_303():
    f = cos(x)**2/(a + b*sin(x)**2)
    F = -x/b + sqrt(a + b)*atan(sqrt(a + b)*tan(x)/sqrt(a))/(sqrt(a)*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_304():
    f = cos(x)/(a + b*sin(x)**2)
    F = atan(sqrt(b)*sin(x)/sqrt(a))/(sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_305():
    f = sec(x)/(a + b*sin(x)**2)
    F = atanh(sin(x))/(a + b) + sqrt(b)*atan(sqrt(b)*sin(x)/sqrt(a))/(sqrt(a)*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_306():
    f = sec(x)**2/(a + b*sin(x)**2)
    F = tan(x)/(a + b) + b*atan(sqrt(a + b)*tan(x)/sqrt(a))/(sqrt(a)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_307():
    f = sec(x)**3/(a + b*sin(x)**2)
    F = tan(x)*sec(x)/(2*a + 2*b) + (a + 3*b)*atanh(sin(x))/(2*(a + b)**2) + b**(sympy.S(3)/2)*atan(sqrt(b)*sin(x)/sqrt(a))/(sqrt(a)*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_308():
    f = sec(x)**4/(a + b*sin(x)**2)
    F = tan(x)**3/(3*a + 3*b) + (a + 2*b)*tan(x)/(a + b)**2 + b**2*atan(sqrt(a + b)*tan(x)/sqrt(a))/(sqrt(a)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_309():
    f = sec(x)**5/(a + b*sin(x)**2)
    F = tan(x)*sec(x)**3/(4*a + 4*b) + (3*a + 7*b)*tan(x)*sec(x)/(8*(a + b)**2) + (3*a**2 + 10*a*b + 15*b**2)*atanh(sin(x))/(8*(a + b)**3) + b**(sympy.S(5)/2)*atan(sqrt(b)*sin(x)/sqrt(a))/(sqrt(a)*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_310():
    f = sec(x)**6/(a + b*sin(x)**2)
    F = tan(x)**5/(5*a + 5*b) + (2*a + 3*b)*tan(x)**3/(3*(a + b)**2) + (a**2 + 3*a*b + 3*b**2)*tan(x)/(a + b)**3 + b**3*atan(sqrt(a + b)*tan(x)/sqrt(a))/(sqrt(a)*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_311():
    f = cos(x)**6/(a + b*sin(x)**2)**2
    F = -sin(x)*cos(x)/(2*b*(a + (a + b)*tan(x)**2)) + x*(4*a + 5*b)/(2*b**3) + (a + b)*(2*a + b)*tan(x)/(2*a*b**2*(a + (a + b)*tan(x)**2)) - (a + b)**(sympy.S(3)/2)*(4*a - b)*atan(sqrt(a + b)*tan(x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_312():
    f = cos(x)**5/(a + b*sin(x)**2)**2
    F = sin(x)/b**2 + (a + b)**2*sin(x)/(2*a*b**2*(a + b*sin(x)**2)) - (a + b)*(3*a - b)*atan(sqrt(b)*sin(x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_313():
    f = cos(x)**4/(a + b*sin(x)**2)**2
    F = x/b**2 + (a + b)*tan(x)/(2*a*b*(a + (a + b)*tan(x)**2)) - sqrt(a + b)*(2*a - b)*atan(sqrt(a + b)*tan(x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_314():
    f = cos(x)**3/(a + b*sin(x)**2)**2
    F = (a + b)*sin(x)/(2*a*b*(a + b*sin(x)**2)) - (a - b)*atan(sqrt(b)*sin(x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_315():
    f = cos(x)**2/(a + b*sin(x)**2)**2
    F = tan(x)/(2*a*(a + (a + b)*tan(x)**2)) + atan(sqrt(a + b)*tan(x)/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_316():
    f = cos(x)/(a + b*sin(x)**2)**2
    F = sin(x)/(2*a*(a + b*sin(x)**2)) + atan(sqrt(b)*sin(x)/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_317():
    f = sec(x)/(a + b*sin(x)**2)**2
    F = atanh(sin(x))/(a + b)**2 + b*sin(x)/(2*a*(a + b)*(a + b*sin(x)**2)) + sqrt(b)*(3*a + b)*atan(sqrt(b)*sin(x)/sqrt(a))/(2*a**(sympy.S(3)/2)*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_318():
    f = sec(x)**2/(a + b*sin(x)**2)**2
    F = tan(x)/(a + b)**2 + b**2*tan(x)/(2*a*(a + b)**2*(a + (a + b)*tan(x)**2)) + b*(4*a + b)*atan(sqrt(a + b)*tan(x)/sqrt(a))/(2*a**(sympy.S(3)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_319():
    f = sec(x)**3/(a + b*sin(x)**2)**2
    F = tan(x)*sec(x)/((a + b*sin(x)**2)*(2*a + 2*b)) + (a + 5*b)*atanh(sin(x))/(2*(a + b)**3) - b*(a - b)*sin(x)/(2*a*(a + b)**2*(a + b*sin(x)**2)) + b**(sympy.S(3)/2)*(5*a + b)*atan(sqrt(b)*sin(x)/sqrt(a))/(2*a**(sympy.S(3)/2)*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_320():
    f = sec(x)**4/(a + b*sin(x)**2)**2
    F = tan(x)**3/(3*(a + b)**2) + (a + 3*b)*tan(x)/(a + b)**3 + b**3*tan(x)/(2*a*(a + b)**3*(a + (a + b)*tan(x)**2)) + b**2*(6*a + b)*atan(sqrt(a + b)*tan(x)/sqrt(a))/(2*a**(sympy.S(3)/2)*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_321():
    f = sqrt(a + b*sin(e + f*x)**2)*cos(e + f*x)**3
    F = a*(a + 4*b)*atanh(sqrt(b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(8*b**(sympy.S(3)/2)*f) + (a + 4*b)*sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)/(8*b*f) - (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)/(4*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_322():
    f = sqrt(a + b*sin(e + f*x)**2)*cos(e + f*x)
    F = a*atanh(sqrt(b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(2*sqrt(b)*f) + sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_323():
    f = sqrt(a + b*sin(e + f*x)**2)*sec(e + f*x)
    F = -sqrt(b)*atanh(sqrt(b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/f + sqrt(a + b)*atanh(sqrt(a + b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_324():
    f = sqrt(a + b*sin(e + f*x)**2)*sec(e + f*x)**3
    F = a*atanh(sqrt(a + b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(2*f*sqrt(a + b)) + sqrt(a + b*sin(e + f*x)**2)*tan(e + f*x)*sec(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_325():
    f = sqrt(a + b*sin(e + f*x)**2)*sec(e + f*x)**5
    F = a*(3*a + 4*b)*atanh(sqrt(a + b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(8*f*(a + b)**(sympy.S(3)/2)) + (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)*sec(e + f*x)**3/(f*(4*a + 4*b)) + sqrt(a + b*sin(e + f*x)**2)*(3*a + 4*b)*tan(e + f*x)*sec(e + f*x)/(f*(8*a + 8*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_326():
    f = sqrt(a + b*sin(e + f*x)**2)
    F = sqrt(a + b*sin(e + f*x)**2)*elliptic_e(e + f*x, -b/a)/(f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_327():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*cos(e + f*x)**3
    F = a**2*(a + 6*b)*atanh(sqrt(b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(16*b**(sympy.S(3)/2)*f) + a*(a + 6*b)*sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)/(16*b*f) + (a + 6*b)*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)/(24*b*f) - (a + b*sin(e + f*x)**2)**(sympy.S(5)/2)*sin(e + f*x)/(6*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_328():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*cos(e + f*x)
    F = 3*a**2*atanh(sqrt(b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(8*sqrt(b)*f) + 3*a*sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)/(8*f) + (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_329():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sec(e + f*x)
    F = -sqrt(b)*(3*a + 2*b)*atanh(sqrt(b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(2*f) - b*sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)/(2*f) + (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_330():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sec(e + f*x)**3
    F = b**(sympy.S(3)/2)*atanh(sqrt(b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/f + (a - 2*b)*sqrt(a + b)*atanh(sqrt(a + b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(2*f) + (a + b)*sqrt(a + b*sin(e + f*x)**2)*tan(e + f*x)*sec(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_331():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sec(e + f*x)**5
    F = 3*a**2*atanh(sqrt(a + b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(8*f*sqrt(a + b)) + 3*a*sqrt(a + b*sin(e + f*x)**2)*tan(e + f*x)*sec(e + f*x)/(8*f) + (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)*sec(e + f*x)**3/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_332():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sec(e + f*x)**7
    F = a**2*(5*a + 6*b)*atanh(sqrt(a + b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(16*f*(a + b)**(sympy.S(3)/2)) + a*sqrt(a + b*sin(e + f*x)**2)*(5*a + 6*b)*tan(e + f*x)*sec(e + f*x)/(f*(16*a + 16*b)) + (a + b*sin(e + f*x)**2)**(sympy.S(5)/2)*tan(e + f*x)*sec(e + f*x)**5/(f*(6*a + 6*b)) + (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*(5*a + 6*b)*tan(e + f*x)*sec(e + f*x)**3/(f*(24*a + 24*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_333():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*cos(e + f*x)**4
    F = a*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*(2*a**2 + 9*a*b - b**2)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(35*b**2*f*sqrt(a + b*sin(e + f*x)**2)) - b*sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)*cos(e + f*x)**5/(7*f) + sqrt(a + b*sin(e + f*x)**2)*(8*a + 2*b)*sin(e + f*x)*cos(e + f*x)**3/(35*f) - sqrt(a + b*sin(e + f*x)**2)*(a**2 - 9*a*b - 2*b**2)*sin(e + f*x)*cos(e + f*x)/(35*b*f) - sqrt(a + b*sin(e + f*x)**2)*(2*a - 2*b)*(a**2 + 6*a*b + b**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(35*b**2*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_334():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*cos(e + f*x)**2
    F = a*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*(3*a - b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(15*b*f*sqrt(a + b*sin(e + f*x)**2)) - b*sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)*cos(e + f*x)**3/(5*f) + sqrt(a + b*sin(e + f*x)**2)*(6*a + 2*b)*sin(e + f*x)*cos(e + f*x)/(15*f) - sqrt(a + b*sin(e + f*x)**2)*(3*a**2 - 7*a*b - 2*b**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(15*b*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_335():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = -a*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*elliptic_f(e + f*x, -b/a)/(3*f*sqrt(a + b*sin(e + f*x)**2)) - b*sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)*cos(e + f*x)/(3*f) + sqrt(a + b*sin(e + f*x)**2)*(4*a + 2*b)*elliptic_e(e + f*x, -b/a)/(3*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_336():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sec(e + f*x)**2
    F = a*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*sqrt(a + b*sin(e + f*x)**2)) + (a + b)*sqrt(a + b*sin(e + f*x)**2)*tan(e + f*x)/f - (a + 2*b)*sqrt(a + b*sin(e + f*x)**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_337():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sec(e + f*x)**4
    F = a*sqrt(1 + b*sin(e + f*x)**2/a)*(2*a - b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(a + b*sin(e + f*x)**2)) + (a + b)*sqrt(a + b*sin(e + f*x)**2)*tan(e + f*x)*sec(e + f*x)**2/(3*f) + sqrt(a + b*sin(e + f*x)**2)*(2*a - 2*b)*tan(e + f*x)/(3*f) - sqrt(a + b*sin(e + f*x)**2)*(2*a - 2*b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_338():
    f = cos(e + f*x)**3/sqrt(a + b*sin(e + f*x)**2)
    F = -sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)/(2*b*f) + (a + 2*b)*atanh(sqrt(b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(2*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_339():
    f = cos(e + f*x)/sqrt(a + b*sin(e + f*x)**2)
    F = atanh(sqrt(b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_340():
    f = sec(e + f*x)/sqrt(a + b*sin(e + f*x)**2)
    F = atanh(sqrt(a + b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_341():
    f = sec(e + f*x)**3/sqrt(a + b*sin(e + f*x)**2)
    F = sqrt(a + b*sin(e + f*x)**2)*tan(e + f*x)*sec(e + f*x)/(f*(2*a + 2*b)) + (a + 2*b)*atanh(sqrt(a + b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(2*f*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_342():
    f = 1/sqrt(a + b*sin(e + f*x)**2)
    F = sqrt(1 + b*sin(e + f*x)**2/a)*elliptic_f(e + f*x, -b/a)/(f*sqrt(a + b*sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_343():
    f = cos(e + f*x)**3/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = -atanh(sqrt(b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(b**(sympy.S(3)/2)*f) + (a + b)*sin(e + f*x)/(a*b*f*sqrt(a + b*sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_344():
    f = cos(e + f*x)/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = sin(e + f*x)/(a*f*sqrt(a + b*sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_345():
    f = sec(e + f*x)/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = atanh(sqrt(a + b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(f*(a + b)**(sympy.S(3)/2)) + b*sin(e + f*x)/(a*f*(a + b)*sqrt(a + b*sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_346():
    f = sec(e + f*x)**3/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = tan(e + f*x)*sec(e + f*x)/(f*sqrt(a + b*sin(e + f*x)**2)*(2*a + 2*b)) + (a + 4*b)*atanh(sqrt(a + b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(2*f*(a + b)**(sympy.S(5)/2)) - b*(a - 2*b)*sin(e + f*x)/(2*a*f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_347():
    f = cos(e + f*x)**6/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = -sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*(8*a + 9*b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*b**3*f*sqrt(a + b*sin(e + f*x)**2)) + (a + b)*sin(e + f*x)*cos(e + f*x)**3/(a*b*f*sqrt(a + b*sin(e + f*x)**2)) + sqrt(a + b*sin(e + f*x)**2)*(4*a + 3*b)*sin(e + f*x)*cos(e + f*x)/(3*a*b**2*f) + sqrt(a + b*sin(e + f*x)**2)*(8*a**2 + 13*a*b + 3*b**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a*b**3*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_348():
    f = cos(e + f*x)**4/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = -sqrt(1 + b*sin(e + f*x)**2/a)*(2*a + 2*b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(b**2*f*sqrt(a + b*sin(e + f*x)**2)) + (a + b)*sin(e + f*x)*cos(e + f*x)/(a*b*f*sqrt(a + b*sin(e + f*x)**2)) + sqrt(a + b*sin(e + f*x)**2)*(2*a + b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(a*b**2*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_349():
    f = cos(e + f*x)**2/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = -sqrt(1 + b*sin(e + f*x)**2/a)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(b*f*sqrt(a + b*sin(e + f*x)**2)) + sin(e + f*x)*cos(e + f*x)/(a*f*sqrt(a + b*sin(e + f*x)**2)) + sqrt(a + b*sin(e + f*x)**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(a*b*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_350():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(-3)/2)
    F = b*sin(e + f*x)*cos(e + f*x)/(a*f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) + sqrt(a + b*sin(e + f*x)**2)*elliptic_e(e + f*x, -b/a)/(a*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_351():
    f = sec(e + f*x)**2/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = sqrt(1 + b*sin(e + f*x)**2/a)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) + tan(e + f*x)/(f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) - b*(a - b)*sin(e + f*x)*cos(e + f*x)/(a*f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2)) - (a - b)*sqrt(a + b*sin(e + f*x)**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(a*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_352():
    f = cos(e + f*x)**5/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = atanh(sqrt(b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(b**(sympy.S(5)/2)*f) + (a + b)*sin(e + f*x)*cos(e + f*x)**2/(3*a*b*f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) - (a + b)*(3*a - 2*b)*sin(e + f*x)/(3*a**2*b**2*f*sqrt(a + b*sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_353():
    f = cos(e + f*x)**3/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = sin(e + f*x)*cos(e + f*x)**2/(3*a*f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) + 2*sin(e + f*x)/(3*a**2*f*sqrt(a + b*sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_354():
    f = cos(e + f*x)/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = sin(e + f*x)/(3*a*f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) + 2*sin(e + f*x)/(3*a**2*f*sqrt(a + b*sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_355():
    f = sec(e + f*x)/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = atanh(sqrt(a + b)*sin(e + f*x)/sqrt(a + b*sin(e + f*x)**2))/(f*(a + b)**(sympy.S(5)/2)) + b*sin(e + f*x)/(3*a*f*(a + b)*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) + b*(5*a + 2*b)*sin(e + f*x)/(3*a**2*f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_356():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(-5)/2)
    F = b*sin(e + f*x)*cos(e + f*x)/(3*a*f*(a + b)*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) - sqrt(1 + b*sin(e + f*x)**2/a)*elliptic_f(e + f*x, -b/a)/(3*a*f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) + 2*b*(2*a + b)*sin(e + f*x)*cos(e + f*x)/(3*a**2*f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2)) + sqrt(a + b*sin(e + f*x)**2)*(4*a + 2*b)*elliptic_e(e + f*x, -b/a)/(3*a**2*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_357():
    f = (d*cos(e + f*x))**m*(a + b*sin(e + f*x)**2)**p
    F = d*(d*cos(e + f*x))**(m - 1)*(a + b*sin(e + f*x)**2)**p*(cos(e + f*x)**2)**(sympy.S.Half - m/2)*sin(e + f*x)*appellf1(sympy.S.Half, -p, sympy.S.Half - m/2, sympy.S(3)/2, -b*sin(e + f*x)**2/a, sin(e + f*x)**2)/(f*(1 + b*sin(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_358():
    f = (a + b*sin(e + f*x)**2)**p*cos(e + f*x)**5
    F = -(a + b*sin(e + f*x)**2)**(p + 1)*sin(e + f*x)*cos(e + f*x)**2/(b*f*(2*p + 5)) - (a + b*sin(e + f*x)**2)**(p + 1)*(3*a + b*(2*p + 7))*sin(e + f*x)/(b**2*f*(2*p + 3)*(2*p + 5)) + (a + b*sin(e + f*x)**2)**p*(3*a**2 + 2*a*b*(2*p + 5) + b**2*(4*p**2 + 16*p + 15))*sin(e + f*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*sin(e + f*x)**2/a)/(b**2*f*(1 + b*sin(e + f*x)**2/a)**p*(2*p + 3)*(2*p + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_359():
    f = (a + b*sin(e + f*x)**2)**p*cos(e + f*x)**3
    F = -(a + b*sin(e + f*x)**2)**(p + 1)*sin(e + f*x)/(b*f*(2*p + 3)) + (a + b*(2*p + 3))*(a + b*sin(e + f*x)**2)**p*sin(e + f*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*sin(e + f*x)**2/a)/(b*f*(1 + b*sin(e + f*x)**2/a)**p*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_360():
    f = (a + b*sin(e + f*x)**2)**p*cos(e + f*x)
    F = (a + b*sin(e + f*x)**2)**p*sin(e + f*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*sin(e + f*x)**2/a)/(f*(1 + b*sin(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_361():
    f = (a + b*sin(e + f*x)**2)**p*sec(e + f*x)
    F = (a + b*sin(e + f*x)**2)**p*sin(e + f*x)*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, sin(e + f*x)**2, -b*sin(e + f*x)**2/a)/(f*(1 + b*sin(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_362():
    f = (a + b*sin(e + f*x)**2)**p*sec(e + f*x)**3
    F = (a + b*sin(e + f*x)**2)**p*sin(e + f*x)*appellf1(sympy.S.Half, 2, -p, sympy.S(3)/2, sin(e + f*x)**2, -b*sin(e + f*x)**2/a)/(f*(1 + b*sin(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_363():
    f = (a + b*sin(e + f*x)**2)**p*cos(e + f*x)**4
    F = (a + b*sin(e + f*x)**2)**p*sqrt(cos(e + f*x)**2)*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S(-3)/2, -p, sympy.S(3)/2, sin(e + f*x)**2, -b*sin(e + f*x)**2/a)/(f*(1 + b*sin(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_364():
    f = (a + b*sin(e + f*x)**2)**p*cos(e + f*x)**2
    F = (a + b*sin(e + f*x)**2)**p*sqrt(cos(e + f*x)**2)*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S(-1)/2, -p, sympy.S(3)/2, sin(e + f*x)**2, -b*sin(e + f*x)**2/a)/(f*(1 + b*sin(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_365():
    f = (a + b*sin(e + f*x)**2)**p
    F = (a + b*sin(e + f*x)**2)**p*sqrt(cos(e + f*x)**2)*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -p, sympy.S(3)/2, sin(e + f*x)**2, -b*sin(e + f*x)**2/a)/(f*(1 + b*sin(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_366():
    f = (a + b*sin(e + f*x)**2)**p*sec(e + f*x)**2
    F = (a + b*sin(e + f*x)**2)**p*sqrt(cos(e + f*x)**2)*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S(3)/2, -p, sympy.S(3)/2, sin(e + f*x)**2, -b*sin(e + f*x)**2/a)/(f*(1 + b*sin(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_367():
    f = (a + b*sin(e + f*x)**2)**p*sec(e + f*x)**4
    F = (a + b*sin(e + f*x)**2)**p*sqrt(cos(e + f*x)**2)*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S(5)/2, -p, sympy.S(3)/2, sin(e + f*x)**2, -b*sin(e + f*x)**2/a)/(f*(1 + b*sin(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_368():
    f = cos(c + d*x)**5/(a + b*sin(c + d*x)**3)
    F = -2*log(a + b*sin(c + d*x)**3)/(3*b*d) + sin(c + d*x)**2/(2*b*d) + sqrt(3)*(a**(sympy.S(4)/3) - b**(sympy.S(4)/3))*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(5)/3)*d) + (a**(sympy.S(4)/3) + b**(sympy.S(4)/3))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(2)/3)*b**(sympy.S(5)/3)*d) - (a**(sympy.S(4)/3) + b**(sympy.S(4)/3))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sin(c + d*x) + b**(sympy.S(2)/3)*sin(c + d*x)**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(5)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_369():
    f = cos(c + d*x)**3/(a + b*sin(c + d*x)**3)
    F = -log(a + b*sin(c + d*x)**3)/(3*b*d) + log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d) - log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sin(c + d*x) + b**(sympy.S(2)/3)*sin(c + d*x)**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d) - sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_370():
    f = cos(c + d*x)/(a + b*sin(c + d*x)**3)
    F = log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d) - log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sin(c + d*x) + b**(sympy.S(2)/3)*sin(c + d*x)**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d) - sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_371():
    f = sec(c + d*x)/(a + b*sin(c + d*x)**3)
    F = -b*log(a + b*sin(c + d*x)**3)/(d*(3*a**2 - 3*b**2)) - log(1 - sin(c + d*x))/(d*(2*a + 2*b)) + log(sin(c + d*x) + 1)/(d*(2*a - 2*b)) - sqrt(3)*b**(sympy.S(1)/3)*(a**(sympy.S(4)/3) - b**(sympy.S(4)/3))*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*d*(a**2 - b**2)) - b**(sympy.S(1)/3)*(a**(sympy.S(4)/3) + b**(sympy.S(4)/3))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(2)/3)*d*(a**2 - b**2)) + b**(sympy.S(1)/3)*(a**(sympy.S(4)/3) + b**(sympy.S(4)/3))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sin(c + d*x) + b**(sympy.S(2)/3)*sin(c + d*x)**2)/(6*a**(sympy.S(2)/3)*d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_372():
    f = sec(c + d*x)**3/(a + b*sin(c + d*x)**3)
    F = b*(a**2 + 2*b**2)*log(a + b*sin(c + d*x)**3)/(3*d*(a**2 - b**2)**2) + (a - 4*b)*log(sin(c + d*x) + 1)/(4*d*(a - b)**2) - 1/(d*(4*a - 4*b)*(sin(c + d*x) + 1)) - (a + 4*b)*log(1 - sin(c + d*x))/(4*d*(a + b)**2) + 1/(d*(1 - sin(c + d*x))*(4*a + 4*b)) - sqrt(3)*b**(sympy.S(5)/3)*(-3*a**(sympy.S(4)/3)*b**(sympy.S(2)/3) + 2*a**2 + b**2)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*d*(a**2 - b**2)**2) + b**(sympy.S(5)/3)*(3*a**(sympy.S(4)/3)*b**(sympy.S(2)/3) + 2*a**2 + b**2)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(2)/3)*d*(a**2 - b**2)**2) - b**(sympy.S(5)/3)*(3*a**(sympy.S(4)/3)*b**(sympy.S(2)/3) + 2*a**2 + b**2)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sin(c + d*x) + b**(sympy.S(2)/3)*sin(c + d*x)**2)/(6*a**(sympy.S(2)/3)*d*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_373():
    f = cos(c + d*x)**4/(a + b*sin(c + d*x)**3)
    F = -2*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*atan((-a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*b**(sympy.S(4)/3)*d*sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) - 2*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(3*b**(sympy.S(4)/3)*d*sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))) + 2*a**(sympy.S(2)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b**(sympy.S(4)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) - cos(c + d*x)/(b*d) + 4*atanh((-(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*b**(sympy.S(2)/3)*d*sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) + 4*atanh(((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*b**(sympy.S(2)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) - 4*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) - 2*atan((-1)**(sympy.S(1)/3)*((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) + 2*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))) + 2*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_374():
    f = cos(c + d*x)**2/(a + b*sin(c + d*x)**3)
    F = 2*atanh((-(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*b**(sympy.S(2)/3)*d*sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) + 2*atanh(((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*b**(sympy.S(2)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) - 2*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) - 2*atan((-1)**(sympy.S(1)/3)*((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) + 2*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))) + 2*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_375():
    f = 1/(a + b*sin(c + d*x)**3)
    F = -2*atan((-1)**(sympy.S(1)/3)*((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) + 2*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))) + 2*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_376():
    f = sec(c + d*x)**2/(a + b*sin(c + d*x)**3)
    F = (-a*sin(c + d*x) + b)*sec(c + d*x)/(d*(-a**2 + b**2)) + 2*(-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)*atan((-a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))**(sympy.S(3)/2)) + 2*(-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))**(sympy.S(3)/2)) - 2*b**(sympy.S(2)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*(a**(sympy.S(2)/3) - b**(sympy.S(2)/3))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_377():
    f = sec(c + d*x)**4/(a + b*sin(c + d*x)**3)
    F = -2*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*b**(sympy.S(8)/3)*atan((-a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(d*sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))*(a**2 - b**2)**2) - 2*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*b**(sympy.S(8)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(d*sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))*(a**2 - b**2)**2) + 2*a**(sympy.S(2)/3)*b**(sympy.S(8)/3)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3))*(a**2 - b**2)**2) - 2*b**(sympy.S(4)/3)*(a**2 + 2*b**2)*atanh((-(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*d*(a**2 - b**2)**2*sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) - 2*b**(sympy.S(4)/3)*(a**2 + 2*b**2)*atanh(((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*d*(a**2 - b**2)**2*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) + 2*b**(sympy.S(4)/3)*(a**2 + 2*b**2)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3))*(a**2 - b**2)**2) - (a - 4*b)*cos(c + d*x)/(4*d*(a - b)**2*(sin(c + d*x) + 1)) - cos(c + d*x)/(d*(12*a - 12*b)*(sin(c + d*x) + 1)) - cos(c + d*x)/(d*(12*a - 12*b)*(sin(c + d*x) + 1)**2) + cos(c + d*x)/(d*(1 - sin(c + d*x))*(12*a + 12*b)) + (a + 4*b)*cos(c + d*x)/(4*d*(1 - sin(c + d*x))*(a + b)**2) + cos(c + d*x)/(d*(1 - sin(c + d*x))**2*(12*a + 12*b)) - 2*b**2*(2*a**2 + b**2)*atan((-a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))*(a**2 - b**2)**2) + 2*b**2*(2*a**2 + b**2)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3))*(a**2 - b**2)**2) + 2*b**2*(2*a**2 + b**2)*atan((a**(sympy.S(1)/3)*tan(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) - b**(sympy.S(2)/3))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_378():
    f = cos(c + d*x)**7/(a + b*sin(c + d*x)**3)**2
    F = -sin(c + d*x)/(b**2*d) - (a**2 + 3*a*b*sin(c + d*x) + 3*b**2*sin(c + d*x)**2 - b**2)*sin(c + d*x)/(3*a*b**2*d*(a + b*sin(c + d*x)**3)) + (-6*a**(sympy.S(4)/3)*b**(sympy.S(2)/3) + 4*a**2 + 2*b**2)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*sin(c + d*x))/(9*a**(sympy.S(5)/3)*b**(sympy.S(7)/3)*d) - (-3*a**(sympy.S(4)/3)*b**(sympy.S(2)/3) + 2*a**2 + b**2)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sin(c + d*x) + b**(sympy.S(2)/3)*sin(c + d*x)**2)/(9*a**(sympy.S(5)/3)*b**(sympy.S(7)/3)*d) - sqrt(3)*(6*a**(sympy.S(4)/3)*b**(sympy.S(2)/3) + 4*a**2 + 2*b**2)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(7)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_379():
    f = cos(c + d*x)**5/(a + b*sin(c + d*x)**3)**2
    F = (-a*sin(c + d*x) - 2*b*sin(c + d*x)**2 + b)*sin(c + d*x)/(3*a*b*d*(a + b*sin(c + d*x)**3)) + (a**(sympy.S(4)/3) - b**(sympy.S(4)/3))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sin(c + d*x) + b**(sympy.S(2)/3)*sin(c + d*x)**2)/(9*a**(sympy.S(5)/3)*b**(sympy.S(5)/3)*d) - (2*a**(sympy.S(4)/3) - 2*b**(sympy.S(4)/3))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*sin(c + d*x))/(9*a**(sympy.S(5)/3)*b**(sympy.S(5)/3)*d) - sqrt(3)*(2*a**(sympy.S(4)/3) + 2*b**(sympy.S(4)/3))*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(5)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_380():
    f = cos(c + d*x)**3/(a + b*sin(c + d*x)**3)**2
    F = (a + b*sin(c + d*x))/(3*a*b*d*(a + b*sin(c + d*x)**3)) + 2*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*sin(c + d*x))/(9*a**(sympy.S(5)/3)*b**(sympy.S(1)/3)*d) - log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sin(c + d*x) + b**(sympy.S(2)/3)*sin(c + d*x)**2)/(9*a**(sympy.S(5)/3)*b**(sympy.S(1)/3)*d) - 2*sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(1)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_381():
    f = cos(c + d*x)/(a + b*sin(c + d*x)**3)**2
    F = sin(c + d*x)/(3*a*d*(a + b*sin(c + d*x)**3)) + 2*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*sin(c + d*x))/(9*a**(sympy.S(5)/3)*b**(sympy.S(1)/3)*d) - log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sin(c + d*x) + b**(sympy.S(2)/3)*sin(c + d*x)**2)/(9*a**(sympy.S(5)/3)*b**(sympy.S(1)/3)*d) - 2*sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(1)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_382():
    f = sec(c + d*x)/(a + b*sin(c + d*x)**3)**2
    F = -2*a*b*log(a + b*sin(c + d*x)**3)/(3*d*(a**2 - b**2)**2) - log(1 - sin(c + d*x))/(2*d*(a + b)**2) + log(sin(c + d*x) + 1)/(2*d*(a - b)**2) + b*(a - (-a*sin(c + d*x) + b)*sin(c + d*x))/(3*a*d*(a + b*sin(c + d*x)**3)*(a**2 - b**2)) - sqrt(3)*b**(sympy.S(1)/3)*(-2*a**(sympy.S(2)/3)*b**(sympy.S(4)/3) + a**2 + b**2)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*d*(a**2 - b**2)**2) - b**(sympy.S(1)/3)*(2*a**(sympy.S(2)/3)*b**(sympy.S(4)/3) + a**2 + b**2)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(1)/3)*d*(a**2 - b**2)**2) + b**(sympy.S(1)/3)*(2*a**(sympy.S(2)/3)*b**(sympy.S(4)/3) + a**2 + b**2)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sin(c + d*x) + b**(sympy.S(2)/3)*sin(c + d*x)**2)/(6*a**(sympy.S(1)/3)*d*(a**2 - b**2)**2) - sqrt(3)*b**(sympy.S(1)/3)*(a**(sympy.S(4)/3) - 2*b**(sympy.S(4)/3))*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*d*(a**2 - b**2)) - b**(sympy.S(1)/3)*(a**(sympy.S(4)/3) + 2*b**(sympy.S(4)/3))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*sin(c + d*x))/(9*a**(sympy.S(5)/3)*d*(a**2 - b**2)) + b**(sympy.S(1)/3)*(a**(sympy.S(4)/3) + 2*b**(sympy.S(4)/3))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sin(c + d*x) + b**(sympy.S(2)/3)*sin(c + d*x)**2)/(18*a**(sympy.S(5)/3)*d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_383():
    f = sec(c + d*x)**3/(a + b*sin(c + d*x)**3)**2
    F = 2*a*b*(a**2 + 5*b**2)*log(a + b*sin(c + d*x)**3)/(3*d*(a**2 - b**2)**3) + (a - 7*b)*log(sin(c + d*x) + 1)/(4*d*(a - b)**3) - (a + 7*b)*log(1 - sin(c + d*x))/(4*d*(a + b)**3) - 1/(4*d*(a - b)**2*(sin(c + d*x) + 1)) + 1/(4*d*(1 - sin(c + d*x))*(a + b)**2) - b*(a*(a**2 + 2*b**2) - b*(2*a**2 - 3*a*b*sin(c + d*x) + b**2)*sin(c + d*x))/(3*a*d*(a + b*sin(c + d*x)**3)*(a**2 - b**2)**2) + b**(sympy.S(5)/3)*(4*a**(sympy.S(2)/3)*(a**2 + 2*b**2) + 3*b**(sympy.S(2)/3)*(3*a**2 + b**2))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(1)/3)*d*(a**2 - b**2)**3) - b**(sympy.S(5)/3)*(4*a**(sympy.S(2)/3)*(a**2 + 2*b**2) + 3*b**(sympy.S(2)/3)*(3*a**2 + b**2))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sin(c + d*x) + b**(sympy.S(2)/3)*sin(c + d*x)**2)/(6*a**(sympy.S(1)/3)*d*(a**2 - b**2)**3) - sqrt(3)*b**(sympy.S(5)/3)*(4*a**(sympy.S(8)/3) + 8*a**(sympy.S(2)/3)*b**2 - 9*a**2*b**(sympy.S(2)/3) - 3*b**(sympy.S(8)/3))*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*d*(a**2 - b**2)**3) - sqrt(3)*b**(sympy.S(5)/3)*(-3*a**(sympy.S(4)/3)*b**(sympy.S(2)/3) + 4*a**2 + 2*b**2)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*sin(c + d*x))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*d*(a**2 - b**2)**2) + b**(sympy.S(5)/3)*(3*a**(sympy.S(4)/3)*b**(sympy.S(2)/3) + 4*a**2 + 2*b**2)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*sin(c + d*x))/(9*a**(sympy.S(5)/3)*d*(a**2 - b**2)**2) - b**(sympy.S(5)/3)*(3*a**(sympy.S(4)/3)*b**(sympy.S(2)/3) + 4*a**2 + 2*b**2)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sin(c + d*x) + b**(sympy.S(2)/3)*sin(c + d*x)**2)/(18*a**(sympy.S(5)/3)*d*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_384():
    f = cos(c + d*x)**4/(a + b*sin(c + d*x)**3)**2
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(4)) * (((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(3)))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_385():
    f = cos(c + d*x)**2/(a + b*sin(c + d*x)**3)**2
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * (((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(3)))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_386():
    f = (a + b*sin(c + d*x)**3)**(-2)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(3)))))**(Integer(2)))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_387():
    f = sec(c + d*x)**2/(a + b*sin(c + d*x)**3)**2
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * (((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(3)))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_388():
    f = sec(c + d*x)**4/(a + b*sin(c + d*x)**3)**2
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**(Integer(4)) * (((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(3)))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_389():
    f = cos(c + d*x)**7/(a - b*sin(c + d*x)**4)
    F = sin(c + d*x)**3/(3*b*d) - 3*sin(c + d*x)/(b*d) - (sqrt(a) - sqrt(b))**3*atanh(b**(sympy.S(1)/4)*sin(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(7)/4)*d) + (sqrt(a) + sqrt(b))**3*atan(b**(sympy.S(1)/4)*sin(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(7)/4)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_390():
    f = cos(c + d*x)**5/(a - b*sin(c + d*x)**4)
    F = -sin(c + d*x)/(b*d) + (sqrt(a) + sqrt(b))**2*atan(b**(sympy.S(1)/4)*sin(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)*d) + (-2*sqrt(a)*sqrt(b) + a + b)*atanh(b**(sympy.S(1)/4)*sin(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_391():
    f = cos(c + d*x)**3/(a - b*sin(c + d*x)**4)
    F = -(sqrt(a) - sqrt(b))*atanh(b**(sympy.S(1)/4)*sin(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*d) + (sqrt(a) + sqrt(b))*atan(b**(sympy.S(1)/4)*sin(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_392():
    f = cos(c + d*x)/(a - b*sin(c + d*x)**4)
    F = atan(b**(sympy.S(1)/4)*sin(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*d) + atanh(b**(sympy.S(1)/4)*sin(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_393():
    f = sec(c + d*x)/(a - b*sin(c + d*x)**4)
    F = atanh(sin(c + d*x))/(d*(a - b)) + b**(sympy.S(1)/4)*atan(b**(sympy.S(1)/4)*sin(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*(sqrt(a) + sqrt(b))) - b**(sympy.S(1)/4)*atanh(b**(sympy.S(1)/4)*sin(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_394():
    f = sec(c + d*x)**3/(a - b*sin(c + d*x)**4)
    F = (a - 5*b)*atanh(sin(c + d*x))/(2*d*(a - b)**2) - 1/(d*(4*a - 4*b)*(sin(c + d*x) + 1)) + 1/(d*(1 - sin(c + d*x))*(4*a - 4*b)) + b**(sympy.S(3)/4)*atan(b**(sympy.S(1)/4)*sin(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*(sqrt(a) + sqrt(b))**2) + b**(sympy.S(3)/4)*atanh(b**(sympy.S(1)/4)*sin(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*(sqrt(a) - sqrt(b))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_395():
    f = sec(c + d*x)**5/(a - b*sin(c + d*x)**4)
    F = -1/(d*(16*a - 16*b)*(sin(c + d*x) + 1)**2) - (3*a - 11*b)/(16*d*(a - b)**2*(sin(c + d*x) + 1)) + (3*a**2 - 6*a*b + 35*b**2)*atanh(sin(c + d*x))/(8*d*(a - b)**3) + (3*a - 11*b)/(16*d*(1 - sin(c + d*x))*(a - b)**2) + 1/(d*(1 - sin(c + d*x))**2*(16*a - 16*b)) + b**(sympy.S(5)/4)*atan(b**(sympy.S(1)/4)*sin(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*(sqrt(a) + sqrt(b))**3) - b**(sympy.S(5)/4)*atanh(b**(sympy.S(1)/4)*sin(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*(sqrt(a) - sqrt(b))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_396():
    f = cos(c + d*x)**10/(a - b*sin(c + d*x)**4)
    F = -17*x/(16*b) - sin(c + d*x)*cos(c + d*x)**5/(6*b*d) - 17*sin(c + d*x)*cos(c + d*x)**3/(24*b*d) - 17*sin(c + d*x)*cos(c + d*x)/(16*b*d) - x*(a + 3*b)/(2*b**2) - x*(4*a + 4*b)/b**2 - (a + 3*b)*sin(c + d*x)*cos(c + d*x)/(2*b**2*d) - (sqrt(a) - sqrt(b))**(sympy.S(9)/2)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(5)/2)*d) + (sqrt(a) + sqrt(b))**(sympy.S(9)/2)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_397():
    f = cos(c + d*x)**8/(a - b*sin(c + d*x)**4)
    F = -11*x/(8*b) - sin(c + d*x)*cos(c + d*x)**3/(4*b*d) - 11*sin(c + d*x)*cos(c + d*x)/(8*b*d) - x*(a + 3*b)/b**2 + (sqrt(a) - sqrt(b))**(sympy.S(7)/2)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**2*d) + (sqrt(a) + sqrt(b))**(sympy.S(7)/2)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_398():
    f = cos(c + d*x)**6/(a - b*sin(c + d*x)**4)
    F = -5*x/(2*b) - sin(c + d*x)*cos(c + d*x)/(2*b*d) - (sqrt(a) - sqrt(b))**(sympy.S(5)/2)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/2)*d) + (sqrt(a) + sqrt(b))**(sympy.S(5)/2)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_399():
    f = cos(c + d*x)**4/(a - b*sin(c + d*x)**4)
    F = -x/b + (sqrt(a) - sqrt(b))**(sympy.S(3)/2)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b*d) + (sqrt(a) + sqrt(b))**(sympy.S(3)/2)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_400():
    f = cos(c + d*x)**2/(a - b*sin(c + d*x)**4)
    F = -sqrt(sqrt(a) - sqrt(b))*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*sqrt(b)*d) + sqrt(sqrt(a) + sqrt(b))*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_401():
    f = sec(c + d*x)**2/(a - b*sin(c + d*x)**4)
    F = tan(c + d*x)/(d*(a - b)) + sqrt(b)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) - sqrt(b)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_402():
    f = sec(c + d*x)**4/(a - b*sin(c + d*x)**4)
    F = (a - 3*b)*tan(c + d*x)/(d*(a - b)**2) + tan(c + d*x)**3/(d*(3*a - 3*b)) + b*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) + b*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_403():
    f = sec(c + d*x)**6/(a - b*sin(c + d*x)**4)
    F = tan(c + d*x)**5/(d*(5*a - 5*b)) + (2*a - 4*b)*tan(c + d*x)**3/(3*d*(a - b)**2) + (a**2 - 3*a*b + 6*b**2)*tan(c + d*x)/(d*(a - b)**3) + b**(sympy.S(3)/2)*atan(sqrt(sqrt(a) + sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(7)/2)) - b**(sympy.S(3)/2)*atan(sqrt(sqrt(a) - sqrt(b))*tan(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_404():
    f = (a + b*sin(e + f*x)**4)**p*cos(e + f*x)**m
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('e') + (Symbol('f') * x))))**(Integer(4)))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_405():
    f = (a + b*sin(e + f*x)**4)**p*cos(e + f*x)**5
    F = -2*(a + b*sin(e + f*x)**4)**p*sin(e + f*x)**3*hyper((sympy.S(3)/4, -p), (sympy.S(7)/4,), -b*sin(e + f*x)**4/a)/(3*f*(1 + b*sin(e + f*x)**4/a)**p) + (a + b*sin(e + f*x)**4)**(p + 1)*sin(e + f*x)/(b*f*(4*p + 5)) - (a - b*(4*p + 5))*(a + b*sin(e + f*x)**4)**p*sin(e + f*x)*hyper((sympy.S(1)/4, -p), (sympy.S(5)/4,), -b*sin(e + f*x)**4/a)/(b*f*(1 + b*sin(e + f*x)**4/a)**p*(4*p + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_406():
    f = (a + b*sin(e + f*x)**4)**p*cos(e + f*x)**3
    F = -(a + b*sin(e + f*x)**4)**p*sin(e + f*x)**3*hyper((sympy.S(3)/4, -p), (sympy.S(7)/4,), -b*sin(e + f*x)**4/a)/(3*f*(1 + b*sin(e + f*x)**4/a)**p) + (a + b*sin(e + f*x)**4)**p*sin(e + f*x)*hyper((sympy.S(1)/4, -p), (sympy.S(5)/4,), -b*sin(e + f*x)**4/a)/(f*(1 + b*sin(e + f*x)**4/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_407():
    f = (a + b*sin(e + f*x)**4)**p*cos(e + f*x)
    F = (a + b*sin(e + f*x)**4)**p*sin(e + f*x)*hyper((sympy.S(1)/4, -p), (sympy.S(5)/4,), -b*sin(e + f*x)**4/a)/(f*(1 + b*sin(e + f*x)**4/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_408():
    f = (a + b*sin(e + f*x)**4)**p*sec(e + f*x)
    F = (a + b*sin(e + f*x)**4)**p*sin(e + f*x)**3*appellf1(sympy.S(3)/4, 1, -p, sympy.S(7)/4, sin(e + f*x)**4, -b*sin(e + f*x)**4/a)/(3*f*(1 + b*sin(e + f*x)**4/a)**p) + (a + b*sin(e + f*x)**4)**p*sin(e + f*x)*appellf1(sympy.S(1)/4, 1, -p, sympy.S(5)/4, sin(e + f*x)**4, -b*sin(e + f*x)**4/a)/(f*(1 + b*sin(e + f*x)**4/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_409():
    f = (a + b*sin(e + f*x)**4)**p*sec(e + f*x)**3
    F = (a + b*sin(e + f*x)**4)**p*sin(e + f*x)**5*appellf1(sympy.S(5)/4, 2, -p, sympy.S(9)/4, sin(e + f*x)**4, -b*sin(e + f*x)**4/a)/(5*f*(1 + b*sin(e + f*x)**4/a)**p) + 2*(a + b*sin(e + f*x)**4)**p*sin(e + f*x)**3*appellf1(sympy.S(3)/4, 2, -p, sympy.S(7)/4, sin(e + f*x)**4, -b*sin(e + f*x)**4/a)/(3*f*(1 + b*sin(e + f*x)**4/a)**p) + (a + b*sin(e + f*x)**4)**p*sin(e + f*x)*appellf1(sympy.S(1)/4, 2, -p, sympy.S(5)/4, sin(e + f*x)**4, -b*sin(e + f*x)**4/a)/(f*(1 + b*sin(e + f*x)**4/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_410():
    f = (a + b*sin(e + f*x)**4)**p*cos(e + f*x)**4
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Integer(4)) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('e') + (Symbol('f') * x))))**(Integer(4)))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_411():
    f = (a + b*sin(e + f*x)**4)**p*cos(e + f*x)**2
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('e') + (Symbol('f') * x))))**(Integer(4)))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_412():
    f = (a + b*sin(e + f*x)**4)**p
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('e') + (Symbol('f') * x))))**(Integer(4)))))**(Symbol('p')), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_413():
    f = (a + b*sin(e + f*x)**4)**p*sec(e + f*x)**2
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('e') + (Symbol('f') * x))))**(Integer(4)))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_414():
    f = (a + b*sin(e + f*x)**4)**p*sec(e + f*x)**4
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(4)) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('e') + (Symbol('f') * x))))**(Integer(4)))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_415():
    f = (a + b*sin(e + f*x)**n)**p*cos(e + f*x)**m
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('e') + (Symbol('f') * x))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_416():
    f = (a + b*sin(e + f*x)**n)**p*cos(e + f*x)**5
    F = (a + b*sin(e + f*x)**n)**p*sin(e + f*x)**5*hyper((-p, 5/n), ((n + 5)/n,), -b*sin(e + f*x)**n/a)/(5*f*(1 + b*sin(e + f*x)**n/a)**p) - 2*(a + b*sin(e + f*x)**n)**p*sin(e + f*x)**3*hyper((-p, 3/n), ((n + 3)/n,), -b*sin(e + f*x)**n/a)/(3*f*(1 + b*sin(e + f*x)**n/a)**p) + (a + b*sin(e + f*x)**n)**p*sin(e + f*x)*hyper((1/n, -p), (1 + 1/n,), -b*sin(e + f*x)**n/a)/(f*(1 + b*sin(e + f*x)**n/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_417():
    f = (a + b*sin(e + f*x)**n)**p*cos(e + f*x)**3
    F = -(a + b*sin(e + f*x)**n)**p*sin(e + f*x)**3*hyper((-p, 3/n), ((n + 3)/n,), -b*sin(e + f*x)**n/a)/(3*f*(1 + b*sin(e + f*x)**n/a)**p) + (a + b*sin(e + f*x)**n)**p*sin(e + f*x)*hyper((1/n, -p), (1 + 1/n,), -b*sin(e + f*x)**n/a)/(f*(1 + b*sin(e + f*x)**n/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_418():
    f = (a + b*sin(e + f*x)**n)**p*cos(e + f*x)
    F = (a + b*sin(e + f*x)**n)**p*sin(e + f*x)*hyper((1/n, -p), (1 + 1/n,), -b*sin(e + f*x)**n/a)/(f*(1 + b*sin(e + f*x)**n/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_419():
    f = (a + b*sin(e + f*x)**n)**p*sec(e + f*x)
    F = sympy.Function('Unintegrable')((sympy.sec((Symbol('e') + (Symbol('f') * x))) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('e') + (Symbol('f') * x))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_420():
    f = (a + b*sin(e + f*x)**n)**p*sec(e + f*x)**3
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('e') + (Symbol('f') * x))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_421():
    f = (a + b*sin(e + f*x)**n)**p*cos(e + f*x)**4
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Integer(4)) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('e') + (Symbol('f') * x))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_422():
    f = (a + b*sin(e + f*x)**n)**p*cos(e + f*x)**2
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('e') + (Symbol('f') * x))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_423():
    f = (a + b*sin(e + f*x)**n)**p
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('e') + (Symbol('f') * x))))**(Symbol('n')))))**(Symbol('p')), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_424():
    f = (a + b*sin(e + f*x)**n)**p*sec(e + f*x)**2
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('e') + (Symbol('f') * x))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_425():
    f = (a + b*sin(e + f*x)**n)**p*sec(e + f*x)**4
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(4)) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('e') + (Symbol('f') * x))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_426():
    f = tan(c + d*x)**7/(a + b*sin(c + d*x)**2)
    F = -a**3*log(a + b*sin(c + d*x)**2)/(2*d*(a + b)**4) + a**3*log(cos(c + d*x))/(d*(a + b)**4) + sec(c + d*x)**6/(d*(6*a + 6*b)) - (3*a + 2*b)*sec(c + d*x)**4/(4*d*(a + b)**2) + (3*a**2 + 3*a*b + b**2)*sec(c + d*x)**2/(2*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_427():
    f = tan(c + d*x)**5/(a + b*sin(c + d*x)**2)
    F = a**2*log(a + b*sin(c + d*x)**2)/(2*d*(a + b)**3) - a**2*log(cos(c + d*x))/(d*(a + b)**3) + sec(c + d*x)**4/(d*(4*a + 4*b)) - (2*a + b)*sec(c + d*x)**2/(2*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_428():
    f = tan(c + d*x)**3/(a + b*sin(c + d*x)**2)
    F = -a*log(a + b*sin(c + d*x)**2)/(2*d*(a + b)**2) + a*log(cos(c + d*x))/(d*(a + b)**2) + sec(c + d*x)**2/(d*(2*a + 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_429():
    f = tan(c + d*x)/(a + b*sin(c + d*x)**2)
    F = log(a + b*sin(c + d*x)**2)/(d*(2*a + 2*b)) - log(cos(c + d*x))/(d*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_430():
    f = cot(c + d*x)/(a + b*sin(c + d*x)**2)
    F = -log(a + b*sin(c + d*x)**2)/(2*a*d) + log(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_431():
    f = cot(c + d*x)**3/(a + b*sin(c + d*x)**2)
    F = -csc(c + d*x)**2/(2*a*d) + (a + b)*log(a + b*sin(c + d*x)**2)/(2*a**2*d) - (a + b)*log(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_432():
    f = cot(c + d*x)**5/(a + b*sin(c + d*x)**2)
    F = -csc(c + d*x)**4/(4*a*d) + (2*a + b)*csc(c + d*x)**2/(2*a**2*d) - (a + b)**2*log(a + b*sin(c + d*x)**2)/(2*a**3*d) + (a + b)**2*log(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_433():
    f = cot(c + d*x)**7/(a + b*sin(c + d*x)**2)
    F = -csc(c + d*x)**6/(6*a*d) + (3*a + b)*csc(c + d*x)**4/(4*a**2*d) - (3*a**2 + 3*a*b + b**2)*csc(c + d*x)**2/(2*a**3*d) + (a + b)**3*log(a + b*sin(c + d*x)**2)/(2*a**4*d) - (a + b)**3*log(sin(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_434():
    f = tan(c + d*x)**8/(a + b*sin(c + d*x)**2)
    F = a**(sympy.S(7)/2)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(d*(a + b)**(sympy.S(9)/2)) - a**3*tan(c + d*x)/(d*(a + b)**4) + a**2*tan(c + d*x)**3/(3*d*(a + b)**3) - a*tan(c + d*x)**5/(5*d*(a + b)**2) + tan(c + d*x)**7/(d*(7*a + 7*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_435():
    f = tan(c + d*x)**6/(a + b*sin(c + d*x)**2)
    F = -a**(sympy.S(5)/2)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(d*(a + b)**(sympy.S(7)/2)) + a**2*tan(c + d*x)/(d*(a + b)**3) - a*tan(c + d*x)**3/(3*d*(a + b)**2) + tan(c + d*x)**5/(d*(5*a + 5*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_436():
    f = tan(c + d*x)**4/(a + b*sin(c + d*x)**2)
    F = a**(sympy.S(3)/2)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(d*(a + b)**(sympy.S(5)/2)) - a*tan(c + d*x)/(d*(a + b)**2) + tan(c + d*x)**3/(d*(3*a + 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_437():
    f = tan(c + d*x)**2/(a + b*sin(c + d*x)**2)
    F = -sqrt(a)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(d*(a + b)**(sympy.S(3)/2)) + tan(c + d*x)/(d*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_438():
    f = cot(c + d*x)**2/(a + b*sin(c + d*x)**2)
    F = -cot(c + d*x)/(a*d) - sqrt(a + b)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_439():
    f = cot(c + d*x)**4/(a + b*sin(c + d*x)**2)
    F = -cot(c + d*x)**3/(3*a*d) + (a + b)*cot(c + d*x)/(a**2*d) + (a + b)**(sympy.S(3)/2)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_440():
    f = cot(c + d*x)**6/(a + b*sin(c + d*x)**2)
    F = -cot(c + d*x)**5/(5*a*d) + (a + b)*cot(c + d*x)**3/(3*a**2*d) - (a + b)**2*cot(c + d*x)/(a**3*d) - (a + b)**(sympy.S(5)/2)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_441():
    f = cot(c + d*x)**8/(a + b*sin(c + d*x)**2)
    F = -cot(c + d*x)**7/(7*a*d) + (a + b)*cot(c + d*x)**5/(5*a**2*d) - (a + b)**2*cot(c + d*x)**3/(3*a**3*d) + (a + b)**3*cot(c + d*x)/(a**4*d) + (a + b)**(sympy.S(7)/2)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(a))/(a**(sympy.S(9)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_442():
    f = sqrt(-a*sin(e + f*x)**2 + a)*tan(e + f*x)**5
    F = a**2/(3*f*(a*cos(e + f*x)**2)**(sympy.S(3)/2)) - 2*a/(f*sqrt(a*cos(e + f*x)**2)) - sqrt(a*cos(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_443():
    f = sqrt(-a*sin(e + f*x)**2 + a)*tan(e + f*x)**3
    F = a/(f*sqrt(a*cos(e + f*x)**2)) + sqrt(a*cos(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_444():
    f = sqrt(-a*sin(e + f*x)**2 + a)*tan(e + f*x)
    F = -sqrt(a*cos(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_445():
    f = sqrt(-a*sin(e + f*x)**2 + a)*cot(e + f*x)
    F = -sqrt(a)*atanh(sqrt(a*cos(e + f*x)**2)/sqrt(a))/f + sqrt(a*cos(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_446():
    f = sqrt(-a*sin(e + f*x)**2 + a)*cot(e + f*x)**3
    F = 3*sqrt(a)*atanh(sqrt(a*cos(e + f*x)**2)/sqrt(a))/(2*f) - 3*sqrt(a*cos(e + f*x)**2)/(2*f) - (a*cos(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**2/(2*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_447():
    f = sqrt(-a*sin(e + f*x)**2 + a)*tan(e + f*x)**6
    F = sqrt(a*cos(e + f*x)**2)*tan(e + f*x)**5/(4*f) - 5*sqrt(a*cos(e + f*x)**2)*tan(e + f*x)**3/(8*f) - 15*sqrt(a*cos(e + f*x)**2)*tan(e + f*x)/(8*f) + 15*sqrt(a*cos(e + f*x)**2)*atanh(sin(e + f*x))*sec(e + f*x)/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_448():
    f = sqrt(-a*sin(e + f*x)**2 + a)*tan(e + f*x)**4
    F = sqrt(a*cos(e + f*x)**2)*tan(e + f*x)**3/(2*f) + 3*sqrt(a*cos(e + f*x)**2)*tan(e + f*x)/(2*f) - 3*sqrt(a*cos(e + f*x)**2)*atanh(sin(e + f*x))*sec(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_449():
    f = sqrt(-a*sin(e + f*x)**2 + a)*tan(e + f*x)**2
    F = -sqrt(a*cos(e + f*x)**2)*tan(e + f*x)/f + sqrt(a*cos(e + f*x)**2)*atanh(sin(e + f*x))*sec(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_450():
    f = sqrt(-a*sin(e + f*x)**2 + a)*cot(e + f*x)**2
    F = -sqrt(a*cos(e + f*x)**2)*tan(e + f*x)/f - sqrt(a*cos(e + f*x)**2)*csc(e + f*x)*sec(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_451():
    f = sqrt(-a*sin(e + f*x)**2 + a)*cot(e + f*x)**4
    F = sqrt(a*cos(e + f*x)**2)*tan(e + f*x)/f - sqrt(a*cos(e + f*x)**2)*csc(e + f*x)**3*sec(e + f*x)/(3*f) + 2*sqrt(a*cos(e + f*x)**2)*csc(e + f*x)*sec(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_452():
    f = sqrt(-a*sin(e + f*x)**2 + a)*cot(e + f*x)**6
    F = -sqrt(a*cos(e + f*x)**2)*tan(e + f*x)/f - sqrt(a*cos(e + f*x)**2)*csc(e + f*x)**5*sec(e + f*x)/(5*f) + sqrt(a*cos(e + f*x)**2)*csc(e + f*x)**3*sec(e + f*x)/f - 3*sqrt(a*cos(e + f*x)**2)*csc(e + f*x)*sec(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_453():
    f = tan(e + f*x)**5/sqrt(-a*sin(e + f*x)**2 + a)
    F = a**2/(5*f*(a*cos(e + f*x)**2)**(sympy.S(5)/2)) - 2*a/(3*f*(a*cos(e + f*x)**2)**(sympy.S(3)/2)) + 1/(f*sqrt(a*cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_454():
    f = tan(e + f*x)**3/sqrt(-a*sin(e + f*x)**2 + a)
    F = a/(3*f*(a*cos(e + f*x)**2)**(sympy.S(3)/2)) - 1/(f*sqrt(a*cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_455():
    f = tan(e + f*x)/sqrt(-a*sin(e + f*x)**2 + a)
    F = 1/(f*sqrt(a*cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_456():
    f = cot(e + f*x)/sqrt(-a*sin(e + f*x)**2 + a)
    F = -atanh(sqrt(a*cos(e + f*x)**2)/sqrt(a))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_457():
    f = cot(e + f*x)**3/sqrt(-a*sin(e + f*x)**2 + a)
    F = -sqrt(a*cos(e + f*x)**2)*csc(e + f*x)**2/(2*a*f) + atanh(sqrt(a*cos(e + f*x)**2)/sqrt(a))/(2*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_458():
    f = tan(e + f*x)**4/sqrt(-a*sin(e + f*x)**2 + a)
    F = 3*cos(e + f*x)*atanh(sin(e + f*x))/(8*f*sqrt(a*cos(e + f*x)**2)) + tan(e + f*x)**3/(4*f*sqrt(a*cos(e + f*x)**2)) - 3*tan(e + f*x)/(8*f*sqrt(a*cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_459():
    f = tan(e + f*x)**2/sqrt(-a*sin(e + f*x)**2 + a)
    F = -cos(e + f*x)*atanh(sin(e + f*x))/(2*f*sqrt(a*cos(e + f*x)**2)) + tan(e + f*x)/(2*f*sqrt(a*cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_460():
    f = cot(e + f*x)**2/sqrt(-a*sin(e + f*x)**2 + a)
    F = -cot(e + f*x)/(f*sqrt(a*cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_461():
    f = cot(e + f*x)**4/sqrt(-a*sin(e + f*x)**2 + a)
    F = -cot(e + f*x)*csc(e + f*x)**2/(3*f*sqrt(a*cos(e + f*x)**2)) + cot(e + f*x)/(f*sqrt(a*cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_462():
    f = cot(e + f*x)**6/sqrt(-a*sin(e + f*x)**2 + a)
    F = -cot(e + f*x)*csc(e + f*x)**4/(5*f*sqrt(a*cos(e + f*x)**2)) + 2*cot(e + f*x)*csc(e + f*x)**2/(3*f*sqrt(a*cos(e + f*x)**2)) - cot(e + f*x)/(f*sqrt(a*cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_463():
    f = tan(e + f*x)**5/(-a*sin(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = a**2/(7*f*(a*cos(e + f*x)**2)**(sympy.S(7)/2)) - 2*a/(5*f*(a*cos(e + f*x)**2)**(sympy.S(5)/2)) + 1/(3*f*(a*cos(e + f*x)**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_464():
    f = tan(e + f*x)**3/(-a*sin(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = a/(5*f*(a*cos(e + f*x)**2)**(sympy.S(5)/2)) - 1/(3*f*(a*cos(e + f*x)**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_465():
    f = tan(e + f*x)/(-a*sin(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = 1/(3*f*(a*cos(e + f*x)**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_466():
    f = cot(e + f*x)/(-a*sin(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = 1/(a*f*sqrt(a*cos(e + f*x)**2)) - atanh(sqrt(a*cos(e + f*x)**2)/sqrt(a))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_467():
    f = cot(e + f*x)**3/(-a*sin(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = -sqrt(a*cos(e + f*x)**2)*csc(e + f*x)**2/(2*a**2*f) - atanh(sqrt(a*cos(e + f*x)**2)/sqrt(a))/(2*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_468():
    f = tan(e + f*x)**2/(-a*sin(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = -cos(e + f*x)*atanh(sin(e + f*x))/(8*a*f*sqrt(a*cos(e + f*x)**2)) + tan(e + f*x)*sec(e + f*x)**2/(4*a*f*sqrt(a*cos(e + f*x)**2)) - tan(e + f*x)/(8*a*f*sqrt(a*cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_469():
    f = cot(e + f*x)**2/(-a*sin(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = cos(e + f*x)*atanh(sin(e + f*x))/(a*f*sqrt(a*cos(e + f*x)**2)) - cot(e + f*x)/(a*f*sqrt(a*cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_470():
    f = cot(e + f*x)**4/(-a*sin(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = -cot(e + f*x)*csc(e + f*x)**2/(3*a*f*sqrt(a*cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_471():
    f = cot(e + f*x)**6/(-a*sin(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = -cot(e + f*x)*csc(e + f*x)**4/(5*a*f*sqrt(a*cos(e + f*x)**2)) + cot(e + f*x)*csc(e + f*x)**2/(3*a*f*sqrt(a*cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_472():
    f = cot(e + f*x)**8/(-a*sin(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = -cot(e + f*x)*csc(e + f*x)**6/(7*a*f*sqrt(a*cos(e + f*x)**2)) + 2*cot(e + f*x)*csc(e + f*x)**4/(5*a*f*sqrt(a*cos(e + f*x)**2)) - cot(e + f*x)*csc(e + f*x)**2/(3*a*f*sqrt(a*cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_473():
    f = sqrt(a + b*sin(e + f*x)**2)*tan(e + f*x)**5
    F = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sec(e + f*x)**4/(f*(4*a + 4*b)) - (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*(8*a + 7*b)*sec(e + f*x)**2/(8*f*(a + b)**2) - sqrt(a + b*sin(e + f*x)**2)*(8*a**2 + 24*a*b + 15*b**2)/(8*f*(a + b)**2) + (8*a**2 + 24*a*b + 15*b**2)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a + b))/(8*f*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_474():
    f = sqrt(a + b*sin(e + f*x)**2)*tan(e + f*x)**3
    F = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*sec(e + f*x)**2/(f*(2*a + 2*b)) + sqrt(a + b*sin(e + f*x)**2)*(2*a + 3*b)/(f*(2*a + 2*b)) - (2*a + 3*b)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a + b))/(2*f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_475():
    f = sqrt(a + b*sin(e + f*x)**2)*tan(e + f*x)
    F = sqrt(a + b)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a + b))/f - sqrt(a + b*sin(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_476():
    f = sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)
    F = -sqrt(a)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a))/f + sqrt(a + b*sin(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_477():
    f = sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)**3
    F = -(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**2/(2*a*f) - sqrt(a + b*sin(e + f*x)**2)*(2*a - b)/(2*a*f) + (2*a - b)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a))/(2*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_478():
    f = sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)**5
    F = -(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**4/(4*a*f) + (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*(8*a + b)*csc(e + f*x)**2/(8*a**2*f) + sqrt(a + b*sin(e + f*x)**2)*(8*a**2 - 8*a*b - b**2)/(8*a**2*f) - (8*a**2 - 8*a*b - b**2)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a))/(8*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_479():
    f = sqrt(a + b*sin(e + f*x)**2)*tan(e + f*x)**4
    F = -4*a*sqrt(1 + b*sin(e + f*x)**2/a)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(a + b*sin(e + f*x)**2)) + sqrt(a + b*sin(e + f*x)**2)*tan(e + f*x)**3/(3*f) - sqrt(a + b*sin(e + f*x)**2)*(3*a + 4*b)*tan(e + f*x)/(f*(3*a + 3*b)) + sqrt(a + b*sin(e + f*x)**2)*(7*a + 8*b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*sqrt(1 + b*sin(e + f*x)**2/a)*(3*a + 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_480():
    f = sqrt(a + b*sin(e + f*x)**2)*tan(e + f*x)**2
    F = a*sqrt(1 + b*sin(e + f*x)**2/a)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*sqrt(a + b*sin(e + f*x)**2)) + sqrt(a + b*sin(e + f*x)**2)*tan(e + f*x)/f - 2*sqrt(a + b*sin(e + f*x)**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_481():
    f = sqrt(a + b*sin(e + f*x)**2)
    F = sqrt(a + b*sin(e + f*x)**2)*elliptic_e(e + f*x, -b/a)/(f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_482():
    f = sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)**2
    F = sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)/f - 2*sqrt(a + b*sin(e + f*x)**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_483():
    f = sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)**4
    F = -sqrt(1 + b*sin(e + f*x)**2/a)*(4*a + 4*b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)**3/(3*f) + sqrt(a + b*sin(e + f*x)**2)*(3*a - b)*cot(e + f*x)/(3*a*f) + sqrt(a + b*sin(e + f*x)**2)*(7*a - b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_484():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)**5
    F = (a + b*sin(e + f*x)**2)**(sympy.S(5)/2)*sec(e + f*x)**4/(f*(4*a + 4*b)) - sqrt(a + b*sin(e + f*x)**2)*(8*a**2 + 40*a*b + 35*b**2)/(f*(8*a + 8*b)) - (a + b*sin(e + f*x)**2)**(sympy.S(5)/2)*(8*a + 9*b)*sec(e + f*x)**2/(8*f*(a + b)**2) - (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*(8*a**2 + 40*a*b + 35*b**2)/(24*f*(a + b)**2) + (8*a**2 + 40*a*b + 35*b**2)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a + b))/(8*f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_485():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)**3
    F = -sqrt(a + b)*(2*a + 5*b)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a + b))/(2*f) + (a + b*sin(e + f*x)**2)**(sympy.S(5)/2)*sec(e + f*x)**2/(f*(2*a + 2*b)) + (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*(2*a + 5*b)/(f*(6*a + 6*b)) + sqrt(a + b*sin(e + f*x)**2)*(2*a + 5*b)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_486():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)
    F = (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a + b))/f - (a + b)*sqrt(a + b*sin(e + f*x)**2)/f - (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_487():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a))/f + a*sqrt(a + b*sin(e + f*x)**2)/f + (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_488():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**3
    F = sqrt(a)*(2*a - 3*b)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a))/(2*f) - sqrt(a + b*sin(e + f*x)**2)*(2*a - 3*b)/(2*f) - (a + b*sin(e + f*x)**2)**(sympy.S(5)/2)*csc(e + f*x)**2/(2*a*f) - (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*(2*a - 3*b)/(6*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_489():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**5
    F = -(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)*csc(e + f*x)**4/(4*a*f) + sqrt(a + b*sin(e + f*x)**2)*(8*a**2 - 24*a*b + 3*b**2)/(8*a*f) + (a + b*sin(e + f*x)**2)**(sympy.S(5)/2)*(8*a - b)*csc(e + f*x)**2/(8*a**2*f) + (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*(8*a**2 - 24*a*b + 3*b**2)/(24*a**2*f) - (8*a**2 - 24*a*b + 3*b**2)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a))/(8*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_490():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)**4
    F = -a*sqrt(1 + b*sin(e + f*x)**2/a)*(5*a + 8*b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(a + b*sin(e + f*x)**2)) - (a + 2*b)*sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)**2*tan(e + f*x)/f + (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)**3/(3*f) - sqrt(a + b*sin(e + f*x)**2)*(3*a + 8*b)*sin(e + f*x)*cos(e + f*x)/(3*f) + sqrt(a + b*sin(e + f*x)**2)*(8*a + 16*b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_491():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)**2
    F = 4*a*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(a + b*sin(e + f*x)**2)) + 4*b*sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)*cos(e + f*x)/(3*f) + (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)/f - sqrt(a + b*sin(e + f*x)**2)*(7*a + 8*b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_492():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = -a*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*elliptic_f(e + f*x, -b/a)/(3*f*sqrt(a + b*sin(e + f*x)**2)) - b*sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)*cos(e + f*x)/(3*f) + sqrt(a + b*sin(e + f*x)**2)*(4*a + 2*b)*elliptic_e(e + f*x, -b/a)/(3*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_493():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**2
    F = 4*a*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(a + b*sin(e + f*x)**2)) + 4*b*sqrt(a + b*sin(e + f*x)**2)*sin(e + f*x)*cos(e + f*x)/(3*f) - (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)/f - sqrt(a + b*sin(e + f*x)**2)*(7*a - b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_494():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**4
    F = -sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*(5*a - 3*b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(a + b*sin(e + f*x)**2)) + (a - b)*sqrt(a + b*sin(e + f*x)**2)*cos(e + f*x)**2*cot(e + f*x)/f - (a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**3/(3*f) + sqrt(a + b*sin(e + f*x)**2)*(3*a - 5*b)*sin(e + f*x)*cos(e + f*x)/(3*f) + sqrt(a + b*sin(e + f*x)**2)*(8*a - 8*b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_495():
    f = tan(e + f*x)**5/sqrt(a + b*sin(e + f*x)**2)
    F = sqrt(a + b*sin(e + f*x)**2)*sec(e + f*x)**4/(f*(4*a + 4*b)) - sqrt(a + b*sin(e + f*x)**2)*(8*a + 5*b)*sec(e + f*x)**2/(8*f*(a + b)**2) + (8*a**2 + 8*a*b + 3*b**2)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a + b))/(8*f*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_496():
    f = tan(e + f*x)**3/sqrt(a + b*sin(e + f*x)**2)
    F = sqrt(a + b*sin(e + f*x)**2)*sec(e + f*x)**2/(f*(2*a + 2*b)) - (2*a + b)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a + b))/(2*f*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_497():
    f = tan(e + f*x)/sqrt(a + b*sin(e + f*x)**2)
    F = atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a + b))/(f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_498():
    f = cot(e + f*x)/sqrt(a + b*sin(e + f*x)**2)
    F = -atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_499():
    f = cot(e + f*x)**3/sqrt(a + b*sin(e + f*x)**2)
    F = -sqrt(a + b*sin(e + f*x)**2)*csc(e + f*x)**2/(2*a*f) + (2*a + b)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a))/(2*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_500():
    f = cot(e + f*x)**5/sqrt(a + b*sin(e + f*x)**2)
    F = -sqrt(a + b*sin(e + f*x)**2)*csc(e + f*x)**4/(4*a*f) + sqrt(a + b*sin(e + f*x)**2)*(8*a + 3*b)*csc(e + f*x)**2/(8*a**2*f) - (8*a**2 + 8*a*b + 3*b**2)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a))/(8*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_501():
    f = tan(e + f*x)**4/sqrt(a + b*sin(e + f*x)**2)
    F = -a*sqrt(1 + b*sin(e + f*x)**2/a)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*sqrt(a + b*sin(e + f*x)**2)*(3*a + 3*b)) + sqrt(a + b*sin(e + f*x)**2)*tan(e + f*x)*sec(e + f*x)**2/(f*(3*a + 3*b)) - sqrt(a + b*sin(e + f*x)**2)*(4*a + 2*b)*tan(e + f*x)/(3*f*(a + b)**2) + sqrt(a + b*sin(e + f*x)**2)*(4*a + 2*b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_502():
    f = tan(e + f*x)**2/sqrt(a + b*sin(e + f*x)**2)
    F = sqrt(a + b*sin(e + f*x)**2)*tan(e + f*x)/(f*(a + b)) - sqrt(a + b*sin(e + f*x)**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_503():
    f = 1/sqrt(a + b*sin(e + f*x)**2)
    F = sqrt(1 + b*sin(e + f*x)**2/a)*elliptic_f(e + f*x, -b/a)/(f*sqrt(a + b*sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_504():
    f = cot(e + f*x)**2/sqrt(a + b*sin(e + f*x)**2)
    F = -sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)/(a*f) - sqrt(a + b*sin(e + f*x)**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(a*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_505():
    f = cot(e + f*x)**4/sqrt(a + b*sin(e + f*x)**2)
    F = -sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a*f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)*csc(e + f*x)**2/(3*a*f) + sqrt(a + b*sin(e + f*x)**2)*(4*a + 2*b)*cot(e + f*x)/(3*a**2*f) + sqrt(a + b*sin(e + f*x)**2)*(4*a + 2*b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a**2*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_506():
    f = tan(e + f*x)**5/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = sec(e + f*x)**4/(f*sqrt(a + b*sin(e + f*x)**2)*(4*a + 4*b)) - (8*a + 3*b)*sec(e + f*x)**2/(8*f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2)) - (8*a**2 - 8*a*b - b**2)/(8*f*(a + b)**3*sqrt(a + b*sin(e + f*x)**2)) + (8*a**2 - 8*a*b - b**2)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a + b))/(8*f*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_507():
    f = tan(e + f*x)**3/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = sec(e + f*x)**2/(f*sqrt(a + b*sin(e + f*x)**2)*(2*a + 2*b)) + (2*a - b)/(2*f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2)) - (2*a - b)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a + b))/(2*f*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_508():
    f = tan(e + f*x)/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = -1/(f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) + atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a + b))/(f*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_509():
    f = cot(e + f*x)/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = 1/(a*f*sqrt(a + b*sin(e + f*x)**2)) - atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_510():
    f = cot(e + f*x)**3/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = -csc(e + f*x)**2/(2*a*f*sqrt(a + b*sin(e + f*x)**2)) - (2*a + 3*b)/(2*a**2*f*sqrt(a + b*sin(e + f*x)**2)) + (2*a + 3*b)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a))/(2*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_511():
    f = cot(e + f*x)**5/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = -csc(e + f*x)**4/(4*a*f*sqrt(a + b*sin(e + f*x)**2)) + (8*a + 5*b)*csc(e + f*x)**2/(8*a**2*f*sqrt(a + b*sin(e + f*x)**2)) + (8*a**2 + 24*a*b + 15*b**2)/(8*a**3*f*sqrt(a + b*sin(e + f*x)**2)) - (8*a**2 + 24*a*b + 15*b**2)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a))/(8*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_512():
    f = tan(e + f*x)**4/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = -4*a*sqrt(1 + b*sin(e + f*x)**2/a)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2)) - 4*a*tan(e + f*x)/(3*f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2)) + b*(7*a - b)*sin(e + f*x)*cos(e + f*x)/(3*f*(a + b)**3*sqrt(a + b*sin(e + f*x)**2)) + tan(e + f*x)*sec(e + f*x)**2/(f*sqrt(a + b*sin(e + f*x)**2)*(3*a + 3*b)) + sqrt(a + b*sin(e + f*x)**2)*(7*a - b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_513():
    f = tan(e + f*x)**2/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = -2*b*sin(e + f*x)*cos(e + f*x)/(f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2)) + sqrt(1 + b*sin(e + f*x)**2/a)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) + tan(e + f*x)/(f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) - 2*sqrt(a + b*sin(e + f*x)**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_514():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(-3)/2)
    F = b*sin(e + f*x)*cos(e + f*x)/(a*f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) + sqrt(a + b*sin(e + f*x)**2)*elliptic_e(e + f*x, -b/a)/(a*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_515():
    f = cot(e + f*x)**2/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = sqrt(1 + b*sin(e + f*x)**2/a)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(a*f*sqrt(a + b*sin(e + f*x)**2)) + cot(e + f*x)/(a*f*sqrt(a + b*sin(e + f*x)**2)) - 2*sqrt(a + b*sin(e + f*x)**2)*cot(e + f*x)/(a**2*f) - 2*sqrt(a + b*sin(e + f*x)**2)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(a**2*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_516():
    f = cot(e + f*x)**4/(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)
    F = (a + b)*cot(e + f*x)*csc(e + f*x)**2/(a*b*f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(1 + b*sin(e + f*x)**2/a)*(4*a + 4*b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a**2*f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*(3*a + 4*b)*cot(e + f*x)*csc(e + f*x)**2/(3*a**2*b*f) + sqrt(a + b*sin(e + f*x)**2)*(7*a + 8*b)*cot(e + f*x)/(3*a**3*f) + sqrt(a + b*sin(e + f*x)**2)*(7*a + 8*b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a**3*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_517():
    f = tan(e + f*x)**5/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = sec(e + f*x)**4/(f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*(4*a + 4*b)) - (8*a + b)*sec(e + f*x)**2/(8*f*(a + b)**2*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) - (8*a**2 - 24*a*b + 3*b**2)/(24*f*(a + b)**3*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) - (8*a**2 - 24*a*b + 3*b**2)/(8*f*(a + b)**4*sqrt(a + b*sin(e + f*x)**2)) + (8*a**2 - 24*a*b + 3*b**2)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a + b))/(8*f*(a + b)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_518():
    f = tan(e + f*x)**3/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = sec(e + f*x)**2/(f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*(2*a + 2*b)) + (2*a - 3*b)/(6*f*(a + b)**2*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) + (2*a - 3*b)/(2*f*(a + b)**3*sqrt(a + b*sin(e + f*x)**2)) - (2*a - 3*b)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a + b))/(2*f*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_519():
    f = tan(e + f*x)/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = -1/(f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*(3*a + 3*b)) - 1/(f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2)) + atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a + b))/(f*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_520():
    f = cot(e + f*x)/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = 1/(3*a*f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) + 1/(a**2*f*sqrt(a + b*sin(e + f*x)**2)) - atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_521():
    f = cot(e + f*x)**3/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = -csc(e + f*x)**2/(2*a*f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) - (2*a + 5*b)/(6*a**2*f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) - (2*a + 5*b)/(2*a**3*f*sqrt(a + b*sin(e + f*x)**2)) + (2*a + 5*b)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a))/(2*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_522():
    f = cot(e + f*x)**5/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = -csc(e + f*x)**4/(4*a*f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) + (8*a + 7*b)*csc(e + f*x)**2/(8*a**2*f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) + (8*a**2 + 40*a*b + 35*b**2)/(24*a**3*f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) + (8*a**2 + 40*a*b + 35*b**2)/(8*a**4*f*sqrt(a + b*sin(e + f*x)**2)) - (8*a**2 + 40*a*b + 35*b**2)*atanh(sqrt(a + b*sin(e + f*x)**2)/sqrt(a))/(8*a**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_523():
    f = tan(e + f*x)**4/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = b*(5*a - 3*b)*sin(e + f*x)*cos(e + f*x)/(3*f*(a + b)**3*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) + b*(8*a - 8*b)*sin(e + f*x)*cos(e + f*x)/(3*f*(a + b)**4*sqrt(a + b*sin(e + f*x)**2)) - sqrt(1 + b*sin(e + f*x)**2/a)*(5*a - 3*b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*(a + b)**3*sqrt(a + b*sin(e + f*x)**2)) + tan(e + f*x)*sec(e + f*x)**2/(f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)*(3*a + 3*b)) - (4*a - 2*b)*tan(e + f*x)/(3*f*(a + b)**2*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) + sqrt(a + b*sin(e + f*x)**2)*(8*a - 8*b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_524():
    f = tan(e + f*x)**2/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = -4*b*sin(e + f*x)*cos(e + f*x)/(3*f*(a + b)**2*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) + 4*sqrt(1 + b*sin(e + f*x)**2/a)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2)) + tan(e + f*x)/(f*(a + b)*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) - b*(7*a - b)*sin(e + f*x)*cos(e + f*x)/(3*a*f*(a + b)**3*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*(7*a - b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_525():
    f = (a + b*sin(e + f*x)**2)**(sympy.S(-5)/2)
    F = b*sin(e + f*x)*cos(e + f*x)/(3*a*f*(a + b)*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) - sqrt(1 + b*sin(e + f*x)**2/a)*elliptic_f(e + f*x, -b/a)/(3*a*f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) + 2*b*(2*a + b)*sin(e + f*x)*cos(e + f*x)/(3*a**2*f*(a + b)**2*sqrt(a + b*sin(e + f*x)**2)) + sqrt(a + b*sin(e + f*x)**2)*(4*a + 2*b)*elliptic_e(e + f*x, -b/a)/(3*a**2*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_526():
    f = cot(e + f*x)**2/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = cot(e + f*x)/(3*a*f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) + 4*sqrt(1 + b*sin(e + f*x)**2/a)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a**2*f*sqrt(a + b*sin(e + f*x)**2)) + (3*a + 4*b)*cot(e + f*x)/(3*a**2*f*(a + b)*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*(7*a + 8*b)*cot(e + f*x)/(3*a**3*f*(a + b)) - sqrt(a + b*sin(e + f*x)**2)*(7*a + 8*b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a**3*f*sqrt(1 + b*sin(e + f*x)**2/a)*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_527():
    f = cot(e + f*x)**4/(a + b*sin(e + f*x)**2)**(sympy.S(5)/2)
    F = (a + b)*cot(e + f*x)*csc(e + f*x)**2/(3*a*b*f*(a + b*sin(e + f*x)**2)**(sympy.S(3)/2)) + (2*a + 6*b)*cot(e + f*x)*csc(e + f*x)**2/(3*a**2*b*f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(1 + b*sin(e + f*x)**2/a)*(5*a + 8*b)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a**3*f*sqrt(a + b*sin(e + f*x)**2)) - sqrt(a + b*sin(e + f*x)**2)*(3*a + 8*b)*cot(e + f*x)*csc(e + f*x)**2/(3*a**3*b*f) + sqrt(a + b*sin(e + f*x)**2)*(8*a + 16*b)*cot(e + f*x)/(3*a**4*f) + sqrt(a + b*sin(e + f*x)**2)*(8*a + 16*b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), -b/a)*sec(e + f*x)/(3*a**4*f*sqrt(1 + b*sin(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_528():
    f = (d*tan(e + f*x))**m*(a + b*sin(e + f*x)**2)**p
    F = (d*tan(e + f*x))**(m + 1)*(a + b*sin(e + f*x)**2)**p*(cos(e + f*x)**2)**(m/2 + sympy.S.Half)*appellf1(m/2 + sympy.S.Half, -p, m/2 + sympy.S.Half, m/2 + sympy.S(3)/2, -b*sin(e + f*x)**2/a, sin(e + f*x)**2)/(d*f*(1 + b*sin(e + f*x)**2/a)**p*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_529():
    f = (a + b*sin(c + d*x)**2)**p*tan(c + d*x)**3
    F = (a + b*sin(c + d*x)**2)**(p + 1)*sec(c + d*x)**2/(d*(2*a + 2*b)) - (a + b*sin(c + d*x)**2)**(p + 1)*(a + b*p + b)*hyper((1, p + 1), (p + 2,), (a + b*sin(c + d*x)**2)/(a + b))/(2*d*(a + b)**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_530():
    f = (a + b*sin(c + d*x)**2)**p*tan(c + d*x)
    F = (a + b*sin(c + d*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), (a + b*sin(c + d*x)**2)/(a + b))/(d*(2*a + 2*b)*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_531():
    f = (a + b*sin(c + d*x)**2)**p*cot(c + d*x)
    F = -(a + b*sin(c + d*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*sin(c + d*x)**2/a)/(2*a*d*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_532():
    f = (a + b*sin(c + d*x)**2)**p*cot(c + d*x)**3
    F = -(a + b*sin(c + d*x)**2)**(p + 1)*csc(c + d*x)**2/(2*a*d) + (a - b*p)*(a + b*sin(c + d*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*sin(c + d*x)**2/a)/(2*a**2*d*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_533():
    f = (a + b*sin(c + d*x)**2)**p*tan(c + d*x)**4
    F = (a + b*sin(c + d*x)**2)**p*sqrt(cos(c + d*x)**2)*sin(c + d*x)**4*tan(c + d*x)*appellf1(sympy.S(5)/2, sympy.S(5)/2, -p, sympy.S(7)/2, sin(c + d*x)**2, -b*sin(c + d*x)**2/a)/(5*d*(1 + b*sin(c + d*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_534():
    f = (a + b*sin(c + d*x)**2)**p*tan(c + d*x)**2
    F = (a + b*sin(c + d*x)**2)**p*sqrt(cos(c + d*x)**2)*sin(c + d*x)**2*tan(c + d*x)*appellf1(sympy.S(3)/2, sympy.S(3)/2, -p, sympy.S(5)/2, sin(c + d*x)**2, -b*sin(c + d*x)**2/a)/(3*d*(1 + b*sin(c + d*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_535():
    f = (a + b*sin(c + d*x)**2)**p*cot(c + d*x)**2
    F = -(a + b*sin(c + d*x)**2)**p*sqrt(cos(c + d*x)**2)*appellf1(sympy.S(-1)/2, sympy.S(-1)/2, -p, sympy.S.Half, sin(c + d*x)**2, -b*sin(c + d*x)**2/a)*csc(c + d*x)*sec(c + d*x)/(d*(1 + b*sin(c + d*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_536():
    f = (a + b*sin(c + d*x)**2)**p*cot(c + d*x)**4
    F = -(a + b*sin(c + d*x)**2)**p*sqrt(cos(c + d*x)**2)*appellf1(sympy.S(-3)/2, sympy.S(-3)/2, -p, sympy.S(-1)/2, sin(c + d*x)**2, -b*sin(c + d*x)**2/a)*csc(c + d*x)**3*sec(c + d*x)/(3*d*(1 + b*sin(c + d*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_537():
    f = cot(x)**3/(a + b*sin(x)**3)
    F = log(a + b*sin(x)**3)/(3*a) - log(sin(x))/a - csc(x)**2/(2*a) - b**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*sin(x))/(3*a**(sympy.S(5)/3)) + b**(sympy.S(2)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sin(x) + b**(sympy.S(2)/3)*sin(x)**2)/(6*a**(sympy.S(5)/3)) + sqrt(3)*b**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*sin(x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_538():
    f = sqrt(a + b*sin(x)**3)*cot(x)
    F = -2*sqrt(a)*atanh(sqrt(a + b*sin(x)**3)/sqrt(a))/3 + 2*sqrt(a + b*sin(x)**3)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_539():
    f = cot(x)/sqrt(a + b*sin(x)**3)
    F = -2*atanh(sqrt(a + b*sin(x)**3)/sqrt(a))/(3*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_540():
    f = sqrt(a + b*sin(c + d*x)**4)*cot(c + d*x)
    F = -sqrt(a)*atanh(sqrt(a + b*sin(c + d*x)**4)/sqrt(a))/(2*d) + sqrt(a + b*sin(c + d*x)**4)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_541():
    f = tan(c + d*x)**3/sqrt(a + b*sin(c + d*x)**4)
    F = -a*atanh((a + b*sin(c + d*x)**2)/(sqrt(a + b)*sqrt(a + b*sin(c + d*x)**4)))/(2*d*(a + b)**(sympy.S(3)/2)) + sqrt(a + b*sin(c + d*x)**4)*sec(c + d*x)**2/(d*(2*a + 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_542():
    f = tan(c + d*x)/sqrt(a + b*sin(c + d*x)**4)
    F = atanh((a + b*sin(c + d*x)**2)/(sqrt(a + b)*sqrt(a + b*sin(c + d*x)**4)))/(2*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_543():
    f = cot(c + d*x)/sqrt(a + b*sin(c + d*x)**4)
    F = -atanh(sqrt(a + b*sin(c + d*x)**4)/sqrt(a))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_544():
    f = cot(c + d*x)**3/sqrt(a + b*sin(c + d*x)**4)
    F = -sqrt(a + b*sin(c + d*x)**4)*csc(c + d*x)**2/(2*a*d) + atanh(sqrt(a + b*sin(c + d*x)**4)/sqrt(a))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_545():
    f = cot(c + d*x)**5/sqrt(a + b*sin(c + d*x)**4)
    F = -sqrt(a + b*sin(c + d*x)**4)*csc(c + d*x)**4/(4*a*d) + sqrt(a + b*sin(c + d*x)**4)*csc(c + d*x)**2/(a*d) - (2*a - b)*atanh(sqrt(a + b*sin(c + d*x)**4)/sqrt(a))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_546():
    f = tan(c + d*x)**2/sqrt(a + b*sin(c + d*x)**4)
    F = -a**(sympy.S(1)/4)*sqrt((2*a*tan(c + d*x)**2 + a + (a + b)*tan(c + d*x)**4)/(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)**2)*(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)*cos(c + d*x)**2*elliptic_e(2*atan((a + b)**(sympy.S(1)/4)*tan(c + d*x)/a**(sympy.S(1)/4)), -sqrt(a)/(2*sqrt(a + b)) + sympy.S.Half)/(d*(a + b)**(sympy.S(3)/4)*sqrt(a + b*sin(c + d*x)**4)) + a**(sympy.S(1)/4)*sqrt((2*a*tan(c + d*x)**2 + a + (a + b)*tan(c + d*x)**4)/(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)**2)*(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)*cos(c + d*x)**2*elliptic_f(2*atan((a + b)**(sympy.S(1)/4)*tan(c + d*x)/a**(sympy.S(1)/4)), -sqrt(a)/(2*sqrt(a + b)) + sympy.S.Half)/(2*d*(a + b)**(sympy.S(3)/4)*sqrt(a + b*sin(c + d*x)**4)) + (2*a*tan(c + d*x)**2 + a + (a + b)*tan(c + d*x)**4)*sin(c + d*x)*cos(c + d*x)/(d*(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)*sqrt(a + b)*sqrt(a + b*sin(c + d*x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_547():
    f = 1/sqrt(a + b*sin(c + d*x)**4)
    F = sqrt((2*a*tan(c + d*x)**2 + a + (a + b)*tan(c + d*x)**4)/(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)**2)*(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)*cos(c + d*x)**2*elliptic_f(2*atan((a + b)**(sympy.S(1)/4)*tan(c + d*x)/a**(sympy.S(1)/4)), -sqrt(a)/(2*sqrt(a + b)) + sympy.S.Half)/(2*a**(sympy.S(1)/4)*d*(a + b)**(sympy.S(1)/4)*sqrt(a + b*sin(c + d*x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_548():
    f = cot(c + d*x)**2/sqrt(a + b*sin(c + d*x)**4)
    F = -(2*a*tan(c + d*x)**2 + a + (a + b)*tan(c + d*x)**4)*cos(c + d*x)**2*cot(c + d*x)/(a*d*sqrt(a + b*sin(c + d*x)**4)) + sqrt(a + b)*(2*a*tan(c + d*x)**2 + a + (a + b)*tan(c + d*x)**4)*sin(c + d*x)*cos(c + d*x)/(a*d*(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)*sqrt(a + b*sin(c + d*x)**4)) - sqrt((2*a*tan(c + d*x)**2 + a + (a + b)*tan(c + d*x)**4)/(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)**2)*(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)*(a + b)**(sympy.S(1)/4)*cos(c + d*x)**2*elliptic_e(2*atan((a + b)**(sympy.S(1)/4)*tan(c + d*x)/a**(sympy.S(1)/4)), -sqrt(a)/(2*sqrt(a + b)) + sympy.S.Half)/(a**(sympy.S(3)/4)*d*sqrt(a + b*sin(c + d*x)**4)) + sqrt((2*a*tan(c + d*x)**2 + a + (a + b)*tan(c + d*x)**4)/(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)**2)*(sqrt(a) + sqrt(a + b)*tan(c + d*x)**2)*(a + b)**(sympy.S(1)/4)*cos(c + d*x)**2*elliptic_f(2*atan((a + b)**(sympy.S(1)/4)*tan(c + d*x)/a**(sympy.S(1)/4)), -sqrt(a)/(2*sqrt(a + b)) + sympy.S.Half)/(2*a**(sympy.S(3)/4)*d*sqrt(a + b*sin(c + d*x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_549():
    f = (a + b*sin(c + d*x)**4)**p*tan(c + d*x)**m
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))**(Symbol('p')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_550():
    f = (a + b*sin(c + d*x)**4)**p*tan(c + d*x)**3
    F = b*(a + b*sin(c + d*x)**4)**p*(2*p + 1)*sin(c + d*x)**2*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*sin(c + d*x)**4/a)/(d*(1 + b*sin(c + d*x)**4/a)**p*(2*a + 2*b)) + (a + b*sin(c + d*x)**4)**(p + 1)*sec(c + d*x)**2/(d*(2*a + 2*b)) - (a + b*sin(c + d*x)**4)**(p + 1)*(a + 2*b*p + b)*hyper((1, p + 1), (p + 2,), (a + b*sin(c + d*x)**4)/(a + b))/(4*d*(a + b)**2*(p + 1)) - (a + b*sin(c + d*x)**4)**p*(a + 2*b*p + b)*sin(c + d*x)**2*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, sin(c + d*x)**4, -b*sin(c + d*x)**4/a)/(d*(1 + b*sin(c + d*x)**4/a)**p*(2*a + 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_551():
    f = (a + b*sin(c + d*x)**4)**p*tan(c + d*x)
    F = (a + b*sin(c + d*x)**4)**(p + 1)*hyper((1, p + 1), (p + 2,), (a + b*sin(c + d*x)**4)/(a + b))/(d*(4*a + 4*b)*(p + 1)) + (a + b*sin(c + d*x)**4)**p*sin(c + d*x)**2*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, sin(c + d*x)**4, -b*sin(c + d*x)**4/a)/(2*d*(1 + b*sin(c + d*x)**4/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_552():
    f = (a + b*sin(c + d*x)**4)**p*cot(c + d*x)
    F = -(a + b*sin(c + d*x)**4)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*sin(c + d*x)**4/a)/(4*a*d*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_553():
    f = (a + b*sin(c + d*x)**4)**p*cot(c + d*x)**3
    F = -(a + b*sin(c + d*x)**4)**p*csc(c + d*x)**2*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*sin(c + d*x)**4/a)/(2*d*(1 + b*sin(c + d*x)**4/a)**p) + (a + b*sin(c + d*x)**4)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*sin(c + d*x)**4/a)/(4*a*d*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_554():
    f = (a + b*sin(c + d*x)**4)**p*tan(c + d*x)**4
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))**(Symbol('p')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_555():
    f = (a + b*sin(c + d*x)**4)**p*tan(c + d*x)**2
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))**(Symbol('p')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_556():
    f = (a + b*sin(c + d*x)**4)**p
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))**(Symbol('p')), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_557():
    f = (a + b*sin(c + d*x)**4)**p*cot(c + d*x)**2
    F = sympy.Function('Unintegrable')(((sympy.cot((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_558():
    f = (a + b*sin(c + d*x)**4)**p*cot(c + d*x)**4
    F = sympy.Function('Unintegrable')(((sympy.cot((Symbol('c') + (Symbol('d') * x))))**(Integer(4)) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(4)))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_559():
    f = (a + b*sin(c + d*x)**n)**3*tan(c + d*x)**m
    F = a**3*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(m + 1)) + 3*a**2*b*(cos(c + d*x)**2)**(m/2 + sympy.S.Half)*sin(c + d*x)**n*tan(c + d*x)**(m + 1)*hyper((m/2 + sympy.S.Half, m/2 + n/2 + sympy.S.Half), (m/2 + n/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*(m + n + 1)) + 3*a*b**2*(cos(c + d*x)**2)**(m/2 + sympy.S.Half)*sin(c + d*x)**(2*n)*tan(c + d*x)**(m + 1)*hyper((m/2 + sympy.S.Half, m/2 + n + sympy.S.Half), (m/2 + n + sympy.S(3)/2,), sin(c + d*x)**2)/(d*(m + 2*n + 1)) + b**3*(cos(c + d*x)**2)**(m/2 + sympy.S.Half)*sin(c + d*x)**(3*n)*tan(c + d*x)**(m + 1)*hyper((m/2 + sympy.S.Half, m/2 + 3*n/2 + sympy.S.Half), (m/2 + 3*n/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*(m + 3*n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_560():
    f = (a + b*sin(c + d*x)**n)**2*tan(c + d*x)**m
    F = a**2*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(m + 1)) + 2*a*b*(cos(c + d*x)**2)**(m/2 + sympy.S.Half)*sin(c + d*x)**n*tan(c + d*x)**(m + 1)*hyper((m/2 + sympy.S.Half, m/2 + n/2 + sympy.S.Half), (m/2 + n/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*(m + n + 1)) + b**2*(cos(c + d*x)**2)**(m/2 + sympy.S.Half)*sin(c + d*x)**(2*n)*tan(c + d*x)**(m + 1)*hyper((m/2 + sympy.S.Half, m/2 + n + sympy.S.Half), (m/2 + n + sympy.S(3)/2,), sin(c + d*x)**2)/(d*(m + 2*n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_561():
    f = (a + b*sin(c + d*x)**n)*tan(c + d*x)**m
    F = a*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(m + 1)) + b*(cos(c + d*x)**2)**(m/2 + sympy.S.Half)*sin(c + d*x)**n*tan(c + d*x)**(m + 1)*hyper((m/2 + sympy.S.Half, m/2 + n/2 + sympy.S.Half), (m/2 + n/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*(m + n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_562():
    f = tan(c + d*x)**m/(a + b*sin(c + d*x)**n)
    F = sympy.Function('Unintegrable')(((sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Symbol('n')))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_563():
    f = tan(c + d*x)**m/(a + b*sin(c + d*x)**n)**2
    F = sympy.Function('Unintegrable')(((sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Symbol('m')) * (((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Symbol('n')))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_564():
    f = sqrt(a + b*sin(x)**n)*cot(x)
    F = -2*sqrt(a)*atanh(sqrt(a + b*sin(x)**n)/sqrt(a))/n + 2*sqrt(a + b*sin(x)**n)/n
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_565():
    f = cot(x)/sqrt(a + b*sin(x)**n)
    F = -2*atanh(sqrt(a + b*sin(x)**n)/sqrt(a))/(sqrt(a)*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_566():
    f = (a + b*sin(c + d*x)**n)**p*tan(c + d*x)**m
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Symbol('n')))))**(Symbol('p')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_567():
    f = (a + b*sin(c + d*x)**n)**p*tan(c + d*x)**3
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Symbol('n')))))**(Symbol('p')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_568():
    f = (a + b*sin(c + d*x)**n)**p*tan(c + d*x)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Symbol('n')))))**(Symbol('p')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_569():
    f = (a + b*sin(c + d*x)**n)**p*cot(c + d*x)
    F = -(a + b*sin(c + d*x)**n)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*sin(c + d*x)**n/a)/(a*d*n*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_570():
    f = (a + b*sin(c + d*x)**n)**p*cot(c + d*x)**3
    F = -(a + b*sin(c + d*x)**n)**p*csc(c + d*x)**2*hyper((-p, -2/n), (-(2 - n)/n,), -b*sin(c + d*x)**n/a)/(2*d*(1 + b*sin(c + d*x)**n/a)**p) + (a + b*sin(c + d*x)**n)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*sin(c + d*x)**n/a)/(a*d*n*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_571():
    f = (a + b*sin(c + d*x)**n)**p*tan(c + d*x)**4
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Symbol('n')))))**(Symbol('p')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_572():
    f = (a + b*sin(c + d*x)**n)**p*tan(c + d*x)**2
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Symbol('n')))))**(Symbol('p')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_573():
    f = (a + b*sin(c + d*x)**n)**p
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Symbol('n')))))**(Symbol('p')), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_574():
    f = (a + b*sin(c + d*x)**n)**p*cot(c + d*x)**2
    F = sympy.Function('Unintegrable')(((sympy.cot((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_575():
    f = (a + b*sin(c + d*x)**n)**p*cot(c + d*x)**4
    F = sympy.Function('Unintegrable')(((sympy.cot((Symbol('c') + (Symbol('d') * x))))**(Integer(4)) * ((Symbol('a') + (Symbol('b') * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_576():
    f = (c*cos(e + f*x))**m*(d*sin(e + f*x))**n*(a + b*sin(e + f*x)**2)**p
    F = c*(c*cos(e + f*x))**(m - 1)*(d*sin(e + f*x))**(n + 1)*(a + b*sin(e + f*x)**2)**p*(cos(e + f*x)**2)**(sympy.S.Half - m/2)*appellf1(n/2 + sympy.S.Half, -p, sympy.S.Half - m/2, n/2 + sympy.S(3)/2, -b*sin(e + f*x)**2/a, sin(e + f*x)**2)/(d*f*(1 + b*sin(e + f*x)**2/a)**p*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_577():
    f = sqrt(a + (b*sin(e + f*x) + c*cos(e + f*x))**2)
    F = sqrt(a + (b*sin(e + f*x) + c*cos(e + f*x))**2)*elliptic_e(e + f*x + atan2(c, b), -(b**2 + c**2)/a)/(f*sqrt(1 + (b*sin(e + f*x) + c*cos(e + f*x))**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_7_d_trig_pow_m_a_plus_b_c_sin_pow_n_pow_p_578():
    f = 1/sqrt(a + (b*sin(e + f*x) + c*cos(e + f*x))**2)
    F = sqrt(1 + (b*sin(e + f*x) + c*cos(e + f*x))**2/a)*elliptic_f(e + f*x + atan2(c, b), -(b**2 + c**2)/a)/(f*sqrt(a + (b*sin(e + f*x) + c*cos(e + f*x))**2))
    assert integrate(f, x) == F

