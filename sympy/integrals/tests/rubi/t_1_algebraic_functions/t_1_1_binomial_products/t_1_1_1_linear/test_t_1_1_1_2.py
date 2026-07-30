"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.1 Binomial products/1.1.1 Linear/1.1.1.2 (a+b x)^m (c+d x)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, m, n, p = symbols('a b c d e m n p')

def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1():
    f = Integer(0)
    F = Integer(0)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_2():
    f = Integer(1)
    F = x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_3():
    f = Integer(5)
    F = 5*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_4():
    f = Integer(-2)
    F = -2*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_5():
    f = sympy.S(-3)/2
    F = -3*x/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_6():
    f = pi
    F = pi*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_7():
    f = a
    F = a*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_8():
    f = 3*a
    F = 3*a*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_9():
    f = pi/sqrt(16 - exp(2))
    F = pi*x/sqrt(16 - exp(2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_10():
    f = x**100
    F = x**101/101
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_11():
    f = x**3
    F = x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_12():
    f = x**2
    F = x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_13():
    f = x
    F = x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_14():
    f = (x)**(Integer(0))
    F = x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_15():
    f = 1/x
    F = log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_16():
    f = x**(-2)
    F = -1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_17():
    f = x**(-3)
    F = -1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_18():
    f = x**(-4)
    F = -1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_19():
    f = x**(-100)
    F = -1/(99*x**99)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_20():
    f = x**(sympy.S(5)/2)
    F = 2*x**(sympy.S(7)/2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_21():
    f = x**(sympy.S(3)/2)
    F = 2*x**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_22():
    f = sqrt(x)
    F = 2*x**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_23():
    f = 1/sqrt(x)
    F = 2*sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_24():
    f = x**(sympy.S(-3)/2)
    F = -2/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_25():
    f = x**(sympy.S(-5)/2)
    F = -2/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_26():
    f = x**(sympy.S(5)/3)
    F = 3*x**(sympy.S(8)/3)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_27():
    f = x**(sympy.S(4)/3)
    F = 3*x**(sympy.S(7)/3)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_28():
    f = x**(sympy.S(2)/3)
    F = 3*x**(sympy.S(5)/3)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_29():
    f = x**(sympy.S(1)/3)
    F = 3*x**(sympy.S(4)/3)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_30():
    f = x**(sympy.S(-1)/3)
    F = 3*x**(sympy.S(2)/3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_31():
    f = x**(sympy.S(-2)/3)
    F = 3*x**(sympy.S(1)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_32():
    f = x**(sympy.S(-4)/3)
    F = -3/x**(sympy.S(1)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_33():
    f = x**(sympy.S(-5)/3)
    F = -3/(2*x**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_34():
    f = x**n
    F = x**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_35():
    f = (b*x)**n
    F = (b*x)**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_36():
    f = 1/(e*(c + d*x) + sqrt(-a))
    F = log(c*e + d*e*x + sqrt(-a))/(d*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_37():
    f = (c + d*(a + b*x))**(sympy.S(5)/2)
    F = 2*(c + d*(a + b*x))**(sympy.S(7)/2)/(7*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_38():
    f = (c + d*(a + b*x))**(sympy.S(3)/2)
    F = 2*(c + d*(a + b*x))**(sympy.S(5)/2)/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_39():
    f = sqrt(c + d*(a + b*x))
    F = 2*(c + d*(a + b*x))**(sympy.S(3)/2)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_40():
    f = 1/sqrt(c + d*(a + b*x))
    F = 2*sqrt(c + d*(a + b*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_41():
    f = (c + d*(a + b*x))**(sympy.S(-3)/2)
    F = -2/(b*d*sqrt(c + d*(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_42():
    f = (c + d*(a + b*x))**(sympy.S(-5)/2)
    F = -2/(3*b*d*(c + d*(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_43():
    f = x**3*(a + b*x)
    F = a*x**4/4 + b*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_44():
    f = x**2*(a + b*x)
    F = a*x**3/3 + b*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_45():
    f = x*(a + b*x)
    F = a*x**2/2 + b*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_46():
    f = a + b*x
    F = a*x + b*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_47():
    f = (a + b*x)/x
    F = a*log(x) + b*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_48():
    f = (a + b*x)/x**2
    F = -a/x + b*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_49():
    f = (a + b*x)/x**3
    F = -(a + b*x)**2/(2*a*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_50():
    f = (a + b*x)/x**4
    F = -a/(3*x**3) - b/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_51():
    f = (a + b*x)/x**5
    F = -a/(4*x**4) - b/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_52():
    f = x**3*(a + b*x)**2
    F = a**2*x**4/4 + 2*a*b*x**5/5 + b**2*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_53():
    f = x**2*(a + b*x)**2
    F = a**2*x**3/3 + a*b*x**4/2 + b**2*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_54():
    f = x*(a + b*x)**2
    F = a**2*x**2/2 + 2*a*b*x**3/3 + b**2*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_55():
    f = (a + b*x)**2
    F = (a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_56():
    f = (a + b*x)**2/x
    F = a**2*log(x) + 2*a*b*x + b**2*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_57():
    f = (a + b*x)**2/x**2
    F = -a**2/x + 2*a*b*log(x) + b**2*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_58():
    f = (a + b*x)**2/x**3
    F = -a**2/(2*x**2) - 2*a*b/x + b**2*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_59():
    f = (a + b*x)**2/x**4
    F = -(a + b*x)**3/(3*a*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_60():
    f = (a + b*x)**2/x**5
    F = -a**2/(4*x**4) - 2*a*b/(3*x**3) - b**2/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_61():
    f = (a + b*x)**2/x**6
    F = -a**2/(5*x**5) - a*b/(2*x**4) - b**2/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_62():
    f = (a + b*x)**2/x**7
    F = -a**2/(6*x**6) - 2*a*b/(5*x**5) - b**2/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_63():
    f = (a + b*x)**2/x**8
    F = -a**2/(7*x**7) - a*b/(3*x**6) - b**2/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_64():
    f = x**4*(a + b*x)**3
    F = a**3*x**5/5 + a**2*b*x**6/2 + 3*a*b**2*x**7/7 + b**3*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_65():
    f = x**3*(a + b*x)**3
    F = a**3*x**4/4 + 3*a**2*b*x**5/5 + a*b**2*x**6/2 + b**3*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_66():
    f = x**2*(a + b*x)**3
    F = a**3*x**3/3 + 3*a**2*b*x**4/4 + 3*a*b**2*x**5/5 + b**3*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_67():
    f = x*(a + b*x)**3
    F = -a*(a + b*x)**4/(4*b**2) + (a + b*x)**5/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_68():
    f = (a + b*x)**3
    F = (a + b*x)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_69():
    f = (a + b*x)**3/x
    F = a**3*log(x) + 3*a**2*b*x + 3*a*b**2*x**2/2 + b**3*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_70():
    f = (a + b*x)**3/x**2
    F = -a**3/x + 3*a**2*b*log(x) + 3*a*b**2*x + b**3*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_71():
    f = (a + b*x)**3/x**3
    F = -a**3/(2*x**2) - 3*a**2*b/x + 3*a*b**2*log(x) + b**3*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_72():
    f = (a + b*x)**3/x**4
    F = -a**3/(3*x**3) - 3*a**2*b/(2*x**2) - 3*a*b**2/x + b**3*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_73():
    f = (a + b*x)**3/x**5
    F = -(a + b*x)**4/(4*a*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_74():
    f = (a + b*x)**3/x**6
    F = -(a + b*x)**4/(5*a*x**5) + b*(a + b*x)**4/(20*a**2*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_75():
    f = (a + b*x)**3/x**7
    F = -a**3/(6*x**6) - 3*a**2*b/(5*x**5) - 3*a*b**2/(4*x**4) - b**3/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_76():
    f = (a + b*x)**3/x**8
    F = -a**3/(7*x**7) - a**2*b/(2*x**6) - 3*a*b**2/(5*x**5) - b**3/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_77():
    f = x**6*(a + b*x)**5
    F = a**5*x**7/7 + 5*a**4*b*x**8/8 + 10*a**3*b**2*x**9/9 + a**2*b**3*x**10 + 5*a*b**4*x**11/11 + b**5*x**12/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_78():
    f = x**5*(a + b*x)**5
    F = a**5*x**6/6 + 5*a**4*b*x**7/7 + 5*a**3*b**2*x**8/4 + 10*a**2*b**3*x**9/9 + a*b**4*x**10/2 + b**5*x**11/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_79():
    f = x**4*(a + b*x)**5
    F = a**5*x**5/5 + 5*a**4*b*x**6/6 + 10*a**3*b**2*x**7/7 + 5*a**2*b**3*x**8/4 + 5*a*b**4*x**9/9 + b**5*x**10/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_80():
    f = x**3*(a + b*x)**5
    F = -a**3*(a + b*x)**6/(6*b**4) + 3*a**2*(a + b*x)**7/(7*b**4) - 3*a*(a + b*x)**8/(8*b**4) + (a + b*x)**9/(9*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_81():
    f = x**2*(a + b*x)**5
    F = a**2*(a + b*x)**6/(6*b**3) - 2*a*(a + b*x)**7/(7*b**3) + (a + b*x)**8/(8*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_82():
    f = x*(a + b*x)**5
    F = -a*(a + b*x)**6/(6*b**2) + (a + b*x)**7/(7*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_83():
    f = (a + b*x)**5
    F = (a + b*x)**6/(6*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_84():
    f = (a + b*x)**5/x
    F = a**5*log(x) + 5*a**4*b*x + 5*a**3*b**2*x**2 + 10*a**2*b**3*x**3/3 + 5*a*b**4*x**4/4 + b**5*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_85():
    f = (a + b*x)**5/x**2
    F = -a**5/x + 5*a**4*b*log(x) + 10*a**3*b**2*x + 5*a**2*b**3*x**2 + 5*a*b**4*x**3/3 + b**5*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_86():
    f = (a + b*x)**5/x**3
    F = -a**5/(2*x**2) - 5*a**4*b/x + 10*a**3*b**2*log(x) + 10*a**2*b**3*x + 5*a*b**4*x**2/2 + b**5*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_87():
    f = (a + b*x)**5/x**4
    F = -a**5/(3*x**3) - 5*a**4*b/(2*x**2) - 10*a**3*b**2/x + 10*a**2*b**3*log(x) + 5*a*b**4*x + b**5*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_88():
    f = (a + b*x)**5/x**5
    F = -a**5/(4*x**4) - 5*a**4*b/(3*x**3) - 5*a**3*b**2/x**2 - 10*a**2*b**3/x + 5*a*b**4*log(x) + b**5*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_89():
    f = (a + b*x)**5/x**6
    F = -a**5/(5*x**5) - 5*a**4*b/(4*x**4) - 10*a**3*b**2/(3*x**3) - 5*a**2*b**3/x**2 - 5*a*b**4/x + b**5*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_90():
    f = (a + b*x)**5/x**7
    F = -(a + b*x)**6/(6*a*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_91():
    f = (a + b*x)**5/x**8
    F = -(a + b*x)**6/(7*a*x**7) + b*(a + b*x)**6/(42*a**2*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_92():
    f = (a + b*x)**5/x**9
    F = -(a + b*x)**6/(8*a*x**8) + b*(a + b*x)**6/(28*a**2*x**7) - b**2*(a + b*x)**6/(168*a**3*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_93():
    f = (a + b*x)**5/x**10
    F = -a**5/(9*x**9) - 5*a**4*b/(8*x**8) - 10*a**3*b**2/(7*x**7) - 5*a**2*b**3/(3*x**6) - a*b**4/x**5 - b**5/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_94():
    f = (a + b*x)**5/x**11
    F = -a**5/(10*x**10) - 5*a**4*b/(9*x**9) - 5*a**3*b**2/(4*x**8) - 10*a**2*b**3/(7*x**7) - 5*a*b**4/(6*x**6) - b**5/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_95():
    f = (a + b*x)**5/x**12
    F = -a**5/(11*x**11) - a**4*b/(2*x**10) - 10*a**3*b**2/(9*x**9) - 5*a**2*b**3/(4*x**8) - 5*a*b**4/(7*x**7) - b**5/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_96():
    f = (a + b*x)**5/x**13
    F = -a**5/(12*x**12) - 5*a**4*b/(11*x**11) - a**3*b**2/x**10 - 10*a**2*b**3/(9*x**9) - 5*a*b**4/(8*x**8) - b**5/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_97():
    f = (a + b*x)**5/x**14
    F = -a**5/(13*x**13) - 5*a**4*b/(12*x**12) - 10*a**3*b**2/(11*x**11) - a**2*b**3/x**10 - 5*a*b**4/(9*x**9) - b**5/(8*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_98():
    f = x**8*(a + b*x)**7
    F = a**7*x**9/9 + 7*a**6*b*x**10/10 + 21*a**5*b**2*x**11/11 + 35*a**4*b**3*x**12/12 + 35*a**3*b**4*x**13/13 + 3*a**2*b**5*x**14/2 + 7*a*b**6*x**15/15 + b**7*x**16/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_99():
    f = x**7*(a + b*x)**7
    F = a**7*x**8/8 + 7*a**6*b*x**9/9 + 21*a**5*b**2*x**10/10 + 35*a**4*b**3*x**11/11 + 35*a**3*b**4*x**12/12 + 21*a**2*b**5*x**13/13 + a*b**6*x**14/2 + b**7*x**15/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_100():
    f = x**6*(a + b*x)**7
    F = a**7*x**7/7 + 7*a**6*b*x**8/8 + 7*a**5*b**2*x**9/3 + 7*a**4*b**3*x**10/2 + 35*a**3*b**4*x**11/11 + 7*a**2*b**5*x**12/4 + 7*a*b**6*x**13/13 + b**7*x**14/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_101():
    f = x**5*(a + b*x)**7
    F = -a**5*(a + b*x)**8/(8*b**6) + 5*a**4*(a + b*x)**9/(9*b**6) - a**3*(a + b*x)**10/b**6 + 10*a**2*(a + b*x)**11/(11*b**6) - 5*a*(a + b*x)**12/(12*b**6) + (a + b*x)**13/(13*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_102():
    f = x**4*(a + b*x)**7
    F = a**4*(a + b*x)**8/(8*b**5) - 4*a**3*(a + b*x)**9/(9*b**5) + 3*a**2*(a + b*x)**10/(5*b**5) - 4*a*(a + b*x)**11/(11*b**5) + (a + b*x)**12/(12*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_103():
    f = x**3*(a + b*x)**7
    F = -a**3*(a + b*x)**8/(8*b**4) + a**2*(a + b*x)**9/(3*b**4) - 3*a*(a + b*x)**10/(10*b**4) + (a + b*x)**11/(11*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_104():
    f = x**2*(a + b*x)**7
    F = a**2*(a + b*x)**8/(8*b**3) - 2*a*(a + b*x)**9/(9*b**3) + (a + b*x)**10/(10*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_105():
    f = x*(a + b*x)**7
    F = -a*(a + b*x)**8/(8*b**2) + (a + b*x)**9/(9*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_106():
    f = (a + b*x)**7
    F = (a + b*x)**8/(8*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_107():
    f = (a + b*x)**7/x
    F = a**7*log(x) + 7*a**6*b*x + 21*a**5*b**2*x**2/2 + 35*a**4*b**3*x**3/3 + 35*a**3*b**4*x**4/4 + 21*a**2*b**5*x**5/5 + 7*a*b**6*x**6/6 + b**7*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_108():
    f = (a + b*x)**7/x**2
    F = -a**7/x + 7*a**6*b*log(x) + 21*a**5*b**2*x + 35*a**4*b**3*x**2/2 + 35*a**3*b**4*x**3/3 + 21*a**2*b**5*x**4/4 + 7*a*b**6*x**5/5 + b**7*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_109():
    f = (a + b*x)**7/x**3
    F = -a**7/(2*x**2) - 7*a**6*b/x + 21*a**5*b**2*log(x) + 35*a**4*b**3*x + 35*a**3*b**4*x**2/2 + 7*a**2*b**5*x**3 + 7*a*b**6*x**4/4 + b**7*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_110():
    f = (a + b*x)**7/x**4
    F = -a**7/(3*x**3) - 7*a**6*b/(2*x**2) - 21*a**5*b**2/x + 35*a**4*b**3*log(x) + 35*a**3*b**4*x + 21*a**2*b**5*x**2/2 + 7*a*b**6*x**3/3 + b**7*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_111():
    f = (a + b*x)**7/x**5
    F = -a**7/(4*x**4) - 7*a**6*b/(3*x**3) - 21*a**5*b**2/(2*x**2) - 35*a**4*b**3/x + 35*a**3*b**4*log(x) + 21*a**2*b**5*x + 7*a*b**6*x**2/2 + b**7*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_112():
    f = (a + b*x)**7/x**6
    F = -a**7/(5*x**5) - 7*a**6*b/(4*x**4) - 7*a**5*b**2/x**3 - 35*a**4*b**3/(2*x**2) - 35*a**3*b**4/x + 21*a**2*b**5*log(x) + 7*a*b**6*x + b**7*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_113():
    f = (a + b*x)**7/x**7
    F = -a**7/(6*x**6) - 7*a**6*b/(5*x**5) - 21*a**5*b**2/(4*x**4) - 35*a**4*b**3/(3*x**3) - 35*a**3*b**4/(2*x**2) - 21*a**2*b**5/x + 7*a*b**6*log(x) + b**7*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_114():
    f = (a + b*x)**7/x**8
    F = -a**7/(7*x**7) - 7*a**6*b/(6*x**6) - 21*a**5*b**2/(5*x**5) - 35*a**4*b**3/(4*x**4) - 35*a**3*b**4/(3*x**3) - 21*a**2*b**5/(2*x**2) - 7*a*b**6/x + b**7*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_115():
    f = (a + b*x)**7/x**9
    F = -(a + b*x)**8/(8*a*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_116():
    f = (a + b*x)**7/x**10
    F = -(a + b*x)**8/(9*a*x**9) + b*(a + b*x)**8/(72*a**2*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_117():
    f = (a + b*x)**7/x**11
    F = -(a + b*x)**8/(10*a*x**10) + b*(a + b*x)**8/(45*a**2*x**9) - b**2*(a + b*x)**8/(360*a**3*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_118():
    f = (a + b*x)**7/x**12
    F = -(a + b*x)**8/(11*a*x**11) + 3*b*(a + b*x)**8/(110*a**2*x**10) - b**2*(a + b*x)**8/(165*a**3*x**9) + b**3*(a + b*x)**8/(1320*a**4*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_119():
    f = (a + b*x)**7/x**13
    F = -(a + b*x)**8/(12*a*x**12) + b*(a + b*x)**8/(33*a**2*x**11) - b**2*(a + b*x)**8/(110*a**3*x**10) + b**3*(a + b*x)**8/(495*a**4*x**9) - b**4*(a + b*x)**8/(3960*a**5*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_120():
    f = (a + b*x)**7/x**14
    F = -a**7/(13*x**13) - 7*a**6*b/(12*x**12) - 21*a**5*b**2/(11*x**11) - 7*a**4*b**3/(2*x**10) - 35*a**3*b**4/(9*x**9) - 21*a**2*b**5/(8*x**8) - a*b**6/x**7 - b**7/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_121():
    f = (a + b*x)**7/x**15
    F = -a**7/(14*x**14) - 7*a**6*b/(13*x**13) - 7*a**5*b**2/(4*x**12) - 35*a**4*b**3/(11*x**11) - 7*a**3*b**4/(2*x**10) - 7*a**2*b**5/(3*x**9) - 7*a*b**6/(8*x**8) - b**7/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_122():
    f = (a + b*x)**7/x**16
    F = -a**7/(15*x**15) - a**6*b/(2*x**14) - 21*a**5*b**2/(13*x**13) - 35*a**4*b**3/(12*x**12) - 35*a**3*b**4/(11*x**11) - 21*a**2*b**5/(10*x**10) - 7*a*b**6/(9*x**9) - b**7/(8*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_123():
    f = x**11*(a + b*x)**10
    F = a**10*x**12/12 + 10*a**9*b*x**13/13 + 45*a**8*b**2*x**14/14 + 8*a**7*b**3*x**15 + 105*a**6*b**4*x**16/8 + 252*a**5*b**5*x**17/17 + 35*a**4*b**6*x**18/3 + 120*a**3*b**7*x**19/19 + 9*a**2*b**8*x**20/4 + 10*a*b**9*x**21/21 + b**10*x**22/22
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_124():
    f = x**10*(a + b*x)**10
    F = a**10*x**11/11 + 5*a**9*b*x**12/6 + 45*a**8*b**2*x**13/13 + 60*a**7*b**3*x**14/7 + 14*a**6*b**4*x**15 + 63*a**5*b**5*x**16/4 + 210*a**4*b**6*x**17/17 + 20*a**3*b**7*x**18/3 + 45*a**2*b**8*x**19/19 + a*b**9*x**20/2 + b**10*x**21/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_125():
    f = x**9*(a + b*x)**10
    F = a**10*x**10/10 + 10*a**9*b*x**11/11 + 15*a**8*b**2*x**12/4 + 120*a**7*b**3*x**13/13 + 15*a**6*b**4*x**14 + 84*a**5*b**5*x**15/5 + 105*a**4*b**6*x**16/8 + 120*a**3*b**7*x**17/17 + 5*a**2*b**8*x**18/2 + 10*a*b**9*x**19/19 + b**10*x**20/20
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_126():
    f = x**8*(a + b*x)**10
    F = a**8*(a + b*x)**11/(11*b**9) - 2*a**7*(a + b*x)**12/(3*b**9) + 28*a**6*(a + b*x)**13/(13*b**9) - 4*a**5*(a + b*x)**14/b**9 + 14*a**4*(a + b*x)**15/(3*b**9) - 7*a**3*(a + b*x)**16/(2*b**9) + 28*a**2*(a + b*x)**17/(17*b**9) - 4*a*(a + b*x)**18/(9*b**9) + (a + b*x)**19/(19*b**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_127():
    f = x**7*(a + b*x)**10
    F = -a**7*(a + b*x)**11/(11*b**8) + 7*a**6*(a + b*x)**12/(12*b**8) - 21*a**5*(a + b*x)**13/(13*b**8) + 5*a**4*(a + b*x)**14/(2*b**8) - 7*a**3*(a + b*x)**15/(3*b**8) + 21*a**2*(a + b*x)**16/(16*b**8) - 7*a*(a + b*x)**17/(17*b**8) + (a + b*x)**18/(18*b**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_128():
    f = x**6*(a + b*x)**10
    F = a**6*(a + b*x)**11/(11*b**7) - a**5*(a + b*x)**12/(2*b**7) + 15*a**4*(a + b*x)**13/(13*b**7) - 10*a**3*(a + b*x)**14/(7*b**7) + a**2*(a + b*x)**15/b**7 - 3*a*(a + b*x)**16/(8*b**7) + (a + b*x)**17/(17*b**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_129():
    f = x**5*(a + b*x)**10
    F = -a**5*(a + b*x)**11/(11*b**6) + 5*a**4*(a + b*x)**12/(12*b**6) - 10*a**3*(a + b*x)**13/(13*b**6) + 5*a**2*(a + b*x)**14/(7*b**6) - a*(a + b*x)**15/(3*b**6) + (a + b*x)**16/(16*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_130():
    f = x**4*(a + b*x)**10
    F = a**4*(a + b*x)**11/(11*b**5) - a**3*(a + b*x)**12/(3*b**5) + 6*a**2*(a + b*x)**13/(13*b**5) - 2*a*(a + b*x)**14/(7*b**5) + (a + b*x)**15/(15*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_131():
    f = x**3*(a + b*x)**10
    F = -a**3*(a + b*x)**11/(11*b**4) + a**2*(a + b*x)**12/(4*b**4) - 3*a*(a + b*x)**13/(13*b**4) + (a + b*x)**14/(14*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_132():
    f = x**2*(a + b*x)**10
    F = a**2*(a + b*x)**11/(11*b**3) - a*(a + b*x)**12/(6*b**3) + (a + b*x)**13/(13*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_133():
    f = x*(a + b*x)**10
    F = -a*(a + b*x)**11/(11*b**2) + (a + b*x)**12/(12*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_134():
    f = (a + b*x)**10
    F = (a + b*x)**11/(11*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_135():
    f = (a + b*x)**10/x
    F = a**10*log(x) + 10*a**9*b*x + 45*a**8*b**2*x**2/2 + 40*a**7*b**3*x**3 + 105*a**6*b**4*x**4/2 + 252*a**5*b**5*x**5/5 + 35*a**4*b**6*x**6 + 120*a**3*b**7*x**7/7 + 45*a**2*b**8*x**8/8 + 10*a*b**9*x**9/9 + b**10*x**10/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_136():
    f = (a + b*x)**10/x**2
    F = -a**10/x + 10*a**9*b*log(x) + 45*a**8*b**2*x + 60*a**7*b**3*x**2 + 70*a**6*b**4*x**3 + 63*a**5*b**5*x**4 + 42*a**4*b**6*x**5 + 20*a**3*b**7*x**6 + 45*a**2*b**8*x**7/7 + 5*a*b**9*x**8/4 + b**10*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_137():
    f = (a + b*x)**10/x**3
    F = -a**10/(2*x**2) - 10*a**9*b/x + 45*a**8*b**2*log(x) + 120*a**7*b**3*x + 105*a**6*b**4*x**2 + 84*a**5*b**5*x**3 + 105*a**4*b**6*x**4/2 + 24*a**3*b**7*x**5 + 15*a**2*b**8*x**6/2 + 10*a*b**9*x**7/7 + b**10*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_138():
    f = (a + b*x)**10/x**4
    F = -a**10/(3*x**3) - 5*a**9*b/x**2 - 45*a**8*b**2/x + 120*a**7*b**3*log(x) + 210*a**6*b**4*x + 126*a**5*b**5*x**2 + 70*a**4*b**6*x**3 + 30*a**3*b**7*x**4 + 9*a**2*b**8*x**5 + 5*a*b**9*x**6/3 + b**10*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_139():
    f = (a + b*x)**10/x**5
    F = -a**10/(4*x**4) - 10*a**9*b/(3*x**3) - 45*a**8*b**2/(2*x**2) - 120*a**7*b**3/x + 210*a**6*b**4*log(x) + 252*a**5*b**5*x + 105*a**4*b**6*x**2 + 40*a**3*b**7*x**3 + 45*a**2*b**8*x**4/4 + 2*a*b**9*x**5 + b**10*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_140():
    f = (a + b*x)**10/x**6
    F = -a**10/(5*x**5) - 5*a**9*b/(2*x**4) - 15*a**8*b**2/x**3 - 60*a**7*b**3/x**2 - 210*a**6*b**4/x + 252*a**5*b**5*log(x) + 210*a**4*b**6*x + 60*a**3*b**7*x**2 + 15*a**2*b**8*x**3 + 5*a*b**9*x**4/2 + b**10*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_141():
    f = (a + b*x)**10/x**7
    F = -a**10/(6*x**6) - 2*a**9*b/x**5 - 45*a**8*b**2/(4*x**4) - 40*a**7*b**3/x**3 - 105*a**6*b**4/x**2 - 252*a**5*b**5/x + 210*a**4*b**6*log(x) + 120*a**3*b**7*x + 45*a**2*b**8*x**2/2 + 10*a*b**9*x**3/3 + b**10*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_142():
    f = (a + b*x)**10/x**8
    F = -a**10/(7*x**7) - 5*a**9*b/(3*x**6) - 9*a**8*b**2/x**5 - 30*a**7*b**3/x**4 - 70*a**6*b**4/x**3 - 126*a**5*b**5/x**2 - 210*a**4*b**6/x + 120*a**3*b**7*log(x) + 45*a**2*b**8*x + 5*a*b**9*x**2 + b**10*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_143():
    f = (a + b*x)**10/x**9
    F = -a**10/(8*x**8) - 10*a**9*b/(7*x**7) - 15*a**8*b**2/(2*x**6) - 24*a**7*b**3/x**5 - 105*a**6*b**4/(2*x**4) - 84*a**5*b**5/x**3 - 105*a**4*b**6/x**2 - 120*a**3*b**7/x + 45*a**2*b**8*log(x) + 10*a*b**9*x + b**10*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_144():
    f = (a + b*x)**10/x**10
    F = -a**10/(9*x**9) - 5*a**9*b/(4*x**8) - 45*a**8*b**2/(7*x**7) - 20*a**7*b**3/x**6 - 42*a**6*b**4/x**5 - 63*a**5*b**5/x**4 - 70*a**4*b**6/x**3 - 60*a**3*b**7/x**2 - 45*a**2*b**8/x + 10*a*b**9*log(x) + b**10*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_145():
    f = (a + b*x)**10/x**11
    F = -a**10/(10*x**10) - 10*a**9*b/(9*x**9) - 45*a**8*b**2/(8*x**8) - 120*a**7*b**3/(7*x**7) - 35*a**6*b**4/x**6 - 252*a**5*b**5/(5*x**5) - 105*a**4*b**6/(2*x**4) - 40*a**3*b**7/x**3 - 45*a**2*b**8/(2*x**2) - 10*a*b**9/x + b**10*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_146():
    f = (a + b*x)**10/x**12
    F = -(a + b*x)**11/(11*a*x**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_147():
    f = (a + b*x)**10/x**13
    F = -(a + b*x)**11/(12*a*x**12) + b*(a + b*x)**11/(132*a**2*x**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_148():
    f = (a + b*x)**10/x**14
    F = -(a + b*x)**11/(13*a*x**13) + b*(a + b*x)**11/(78*a**2*x**12) - b**2*(a + b*x)**11/(858*a**3*x**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_149():
    f = (a + b*x)**10/x**15
    F = -(a + b*x)**11/(14*a*x**14) + 3*b*(a + b*x)**11/(182*a**2*x**13) - b**2*(a + b*x)**11/(364*a**3*x**12) + b**3*(a + b*x)**11/(4004*a**4*x**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_150():
    f = (a + b*x)**10/x**16
    F = -(a + b*x)**11/(15*a*x**15) + 2*b*(a + b*x)**11/(105*a**2*x**14) - 2*b**2*(a + b*x)**11/(455*a**3*x**13) + b**3*(a + b*x)**11/(1365*a**4*x**12) - b**4*(a + b*x)**11/(15015*a**5*x**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_151():
    f = (a + b*x)**10/x**17
    F = -(a + b*x)**11/(16*a*x**16) + b*(a + b*x)**11/(48*a**2*x**15) - b**2*(a + b*x)**11/(168*a**3*x**14) + b**3*(a + b*x)**11/(728*a**4*x**13) - b**4*(a + b*x)**11/(4368*a**5*x**12) + b**5*(a + b*x)**11/(48048*a**6*x**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_152():
    f = (a + b*x)**10/x**18
    F = -(a + b*x)**11/(17*a*x**17) + 3*b*(a + b*x)**11/(136*a**2*x**16) - b**2*(a + b*x)**11/(136*a**3*x**15) + b**3*(a + b*x)**11/(476*a**4*x**14) - 3*b**4*(a + b*x)**11/(6188*a**5*x**13) + b**5*(a + b*x)**11/(12376*a**6*x**12) - b**6*(a + b*x)**11/(136136*a**7*x**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_153():
    f = (a + b*x)**10/x**19
    F = -a**10/(18*x**18) - 10*a**9*b/(17*x**17) - 45*a**8*b**2/(16*x**16) - 8*a**7*b**3/x**15 - 15*a**6*b**4/x**14 - 252*a**5*b**5/(13*x**13) - 35*a**4*b**6/(2*x**12) - 120*a**3*b**7/(11*x**11) - 9*a**2*b**8/(2*x**10) - 10*a*b**9/(9*x**9) - b**10/(8*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_154():
    f = (a + b*x)**10/x**20
    F = -a**10/(19*x**19) - 5*a**9*b/(9*x**18) - 45*a**8*b**2/(17*x**17) - 15*a**7*b**3/(2*x**16) - 14*a**6*b**4/x**15 - 18*a**5*b**5/x**14 - 210*a**4*b**6/(13*x**13) - 10*a**3*b**7/x**12 - 45*a**2*b**8/(11*x**11) - a*b**9/x**10 - b**10/(9*x**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_155():
    f = c*(a + b*x)
    F = c*(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_156():
    f = (a + b*x)*(c + d)/e
    F = (a + b*x)**2*(c + d)/(2*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_157():
    f = x**5/(a + b*x)
    F = -a**5*log(a + b*x)/b**6 + a**4*x/b**5 - a**3*x**2/(2*b**4) + a**2*x**3/(3*b**3) - a*x**4/(4*b**2) + x**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_158():
    f = x**4/(a + b*x)
    F = a**4*log(a + b*x)/b**5 - a**3*x/b**4 + a**2*x**2/(2*b**3) - a*x**3/(3*b**2) + x**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_159():
    f = x**3/(a + b*x)
    F = -a**3*log(a + b*x)/b**4 + a**2*x/b**3 - a*x**2/(2*b**2) + x**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_160():
    f = x**2/(a + b*x)
    F = a**2*log(a + b*x)/b**3 - a*x/b**2 + x**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_161():
    f = x/(a + b*x)
    F = -a*log(a + b*x)/b**2 + x/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_162():
    f = 1/(a + b*x)
    F = log(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_163():
    f = 1/(x*(a + b*x))
    F = log(x)/a - log(a + b*x)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_164():
    f = 1/(x**2*(a + b*x))
    F = -1/(a*x) - b*log(x)/a**2 + b*log(a + b*x)/a**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_165():
    f = 1/(x**3*(a + b*x))
    F = -1/(2*a*x**2) + b/(a**2*x) + b**2*log(x)/a**3 - b**2*log(a + b*x)/a**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_166():
    f = 1/(x**4*(a + b*x))
    F = -1/(3*a*x**3) + b/(2*a**2*x**2) - b**2/(a**3*x) - b**3*log(x)/a**4 + b**3*log(a + b*x)/a**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_167():
    f = 1/(x**5*(a + b*x))
    F = -1/(4*a*x**4) + b/(3*a**2*x**3) - b**2/(2*a**3*x**2) + b**3/(a**4*x) + b**4*log(x)/a**5 - b**4*log(a + b*x)/a**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_168():
    f = x**6/(a + b*x)**2
    F = -a**6/(b**7*(a + b*x)) - 6*a**5*log(a + b*x)/b**7 + 5*a**4*x/b**6 - 2*a**3*x**2/b**5 + a**2*x**3/b**4 - a*x**4/(2*b**3) + x**5/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_169():
    f = x**5/(a + b*x)**2
    F = a**5/(b**6*(a + b*x)) + 5*a**4*log(a + b*x)/b**6 - 4*a**3*x/b**5 + 3*a**2*x**2/(2*b**4) - 2*a*x**3/(3*b**3) + x**4/(4*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_170():
    f = x**4/(a + b*x)**2
    F = -a**4/(b**5*(a + b*x)) - 4*a**3*log(a + b*x)/b**5 + 3*a**2*x/b**4 - a*x**2/b**3 + x**3/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_171():
    f = x**3/(a + b*x)**2
    F = a**3/(b**4*(a + b*x)) + 3*a**2*log(a + b*x)/b**4 - 2*a*x/b**3 + x**2/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_172():
    f = x**2/(a + b*x)**2
    F = -a**2/(b**3*(a + b*x)) - 2*a*log(a + b*x)/b**3 + x/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_173():
    f = x/(a + b*x)**2
    F = a/(b**2*(a + b*x)) + log(a + b*x)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_174():
    f = (a + b*x)**(-2)
    F = -1/(b*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_175():
    f = 1/(x*(a + b*x)**2)
    F = 1/(a*(a + b*x)) + log(x)/a**2 - log(a + b*x)/a**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_176():
    f = 1/(x**2*(a + b*x)**2)
    F = -b/(a**2*(a + b*x)) - 1/(a**2*x) - 2*b*log(x)/a**3 + 2*b*log(a + b*x)/a**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_177():
    f = 1/(x**3*(a + b*x)**2)
    F = -1/(2*a**2*x**2) + b**2/(a**3*(a + b*x)) + 2*b/(a**3*x) + 3*b**2*log(x)/a**4 - 3*b**2*log(a + b*x)/a**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_178():
    f = 1/(x**4*(a + b*x)**2)
    F = -1/(3*a**2*x**3) + b/(a**3*x**2) - b**3/(a**4*(a + b*x)) - 3*b**2/(a**4*x) - 4*b**3*log(x)/a**5 + 4*b**3*log(a + b*x)/a**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_179():
    f = 1/(x**5*(a + b*x)**2)
    F = -1/(4*a**2*x**4) + 2*b/(3*a**3*x**3) - 3*b**2/(2*a**4*x**2) + b**4/(a**5*(a + b*x)) + 4*b**3/(a**5*x) + 5*b**4*log(x)/a**6 - 5*b**4*log(a + b*x)/a**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_180():
    f = x**7/(a + b*x)**3
    F = a**7/(2*b**8*(a + b*x)**2) - 7*a**6/(b**8*(a + b*x)) - 21*a**5*log(a + b*x)/b**8 + 15*a**4*x/b**7 - 5*a**3*x**2/b**6 + 2*a**2*x**3/b**5 - 3*a*x**4/(4*b**4) + x**5/(5*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_181():
    f = x**6/(a + b*x)**3
    F = -a**6/(2*b**7*(a + b*x)**2) + 6*a**5/(b**7*(a + b*x)) + 15*a**4*log(a + b*x)/b**7 - 10*a**3*x/b**6 + 3*a**2*x**2/b**5 - a*x**3/b**4 + x**4/(4*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_182():
    f = x**5/(a + b*x)**3
    F = a**5/(2*b**6*(a + b*x)**2) - 5*a**4/(b**6*(a + b*x)) - 10*a**3*log(a + b*x)/b**6 + 6*a**2*x/b**5 - 3*a*x**2/(2*b**4) + x**3/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_183():
    f = x**4/(a + b*x)**3
    F = -a**4/(2*b**5*(a + b*x)**2) + 4*a**3/(b**5*(a + b*x)) + 6*a**2*log(a + b*x)/b**5 - 3*a*x/b**4 + x**2/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_184():
    f = x**3/(a + b*x)**3
    F = a**3/(2*b**4*(a + b*x)**2) - 3*a**2/(b**4*(a + b*x)) - 3*a*log(a + b*x)/b**4 + x/b**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_185():
    f = x**2/(a + b*x)**3
    F = -a**2/(2*b**3*(a + b*x)**2) + 2*a/(b**3*(a + b*x)) + log(a + b*x)/b**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_186():
    f = x/(a + b*x)**3
    F = x**2/(2*a*(a + b*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_187():
    f = (a + b*x)**(-3)
    F = -1/(2*b*(a + b*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_188():
    f = 1/(x*(a + b*x)**3)
    F = 1/(2*a*(a + b*x)**2) + 1/(a**2*(a + b*x)) + log(x)/a**3 - log(a + b*x)/a**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_189():
    f = 1/(x**2*(a + b*x)**3)
    F = -b/(2*a**2*(a + b*x)**2) - 2*b/(a**3*(a + b*x)) - 1/(a**3*x) - 3*b*log(x)/a**4 + 3*b*log(a + b*x)/a**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_190():
    f = 1/(x**3*(a + b*x)**3)
    F = b**2/(2*a**3*(a + b*x)**2) - 1/(2*a**3*x**2) + 3*b**2/(a**4*(a + b*x)) + 3*b/(a**4*x) + 6*b**2*log(x)/a**5 - 6*b**2*log(a + b*x)/a**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_191():
    f = 1/(x**4*(a + b*x)**3)
    F = -1/(3*a**3*x**3) - b**3/(2*a**4*(a + b*x)**2) + 3*b/(2*a**4*x**2) - 4*b**3/(a**5*(a + b*x)) - 6*b**2/(a**5*x) - 10*b**3*log(x)/a**6 + 10*b**3*log(a + b*x)/a**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_192():
    f = 1/(x**5*(a + b*x)**3)
    F = -1/(4*a**3*x**4) + b/(a**4*x**3) + b**4/(2*a**5*(a + b*x)**2) - 3*b**2/(a**5*x**2) + 5*b**4/(a**6*(a + b*x)) + 10*b**3/(a**6*x) + 15*b**4*log(x)/a**7 - 15*b**4*log(a + b*x)/a**7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_193():
    f = x**8/(a + b*x)**4
    F = -a**8/(3*b**9*(a + b*x)**3) + 4*a**7/(b**9*(a + b*x)**2) - 28*a**6/(b**9*(a + b*x)) - 56*a**5*log(a + b*x)/b**9 + 35*a**4*x/b**8 - 10*a**3*x**2/b**7 + 10*a**2*x**3/(3*b**6) - a*x**4/b**5 + x**5/(5*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_194():
    f = x**7/(a + b*x)**4
    F = a**7/(3*b**8*(a + b*x)**3) - 7*a**6/(2*b**8*(a + b*x)**2) + 21*a**5/(b**8*(a + b*x)) + 35*a**4*log(a + b*x)/b**8 - 20*a**3*x/b**7 + 5*a**2*x**2/b**6 - 4*a*x**3/(3*b**5) + x**4/(4*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_195():
    f = x**6/(a + b*x)**4
    F = -a**6/(3*b**7*(a + b*x)**3) + 3*a**5/(b**7*(a + b*x)**2) - 15*a**4/(b**7*(a + b*x)) - 20*a**3*log(a + b*x)/b**7 + 10*a**2*x/b**6 - 2*a*x**2/b**5 + x**3/(3*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_196():
    f = x**5/(a + b*x)**4
    F = a**5/(3*b**6*(a + b*x)**3) - 5*a**4/(2*b**6*(a + b*x)**2) + 10*a**3/(b**6*(a + b*x)) + 10*a**2*log(a + b*x)/b**6 - 4*a*x/b**5 + x**2/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_197():
    f = x**4/(a + b*x)**4
    F = -a**4/(3*b**5*(a + b*x)**3) + 2*a**3/(b**5*(a + b*x)**2) - 6*a**2/(b**5*(a + b*x)) - 4*a*log(a + b*x)/b**5 + x/b**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_198():
    f = x**3/(a + b*x)**4
    F = a**3/(3*b**4*(a + b*x)**3) - 3*a**2/(2*b**4*(a + b*x)**2) + 3*a/(b**4*(a + b*x)) + log(a + b*x)/b**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_199():
    f = x**2/(a + b*x)**4
    F = x**3/(3*a*(a + b*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_200():
    f = x/(a + b*x)**4
    F = a/(3*b**2*(a + b*x)**3) - 1/(2*b**2*(a + b*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_201():
    f = (a + b*x)**(-4)
    F = -1/(3*b*(a + b*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_202():
    f = 1/(x*(a + b*x)**4)
    F = 1/(3*a*(a + b*x)**3) + 1/(2*a**2*(a + b*x)**2) + 1/(a**3*(a + b*x)) + log(x)/a**4 - log(a + b*x)/a**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_203():
    f = 1/(x**2*(a + b*x)**4)
    F = -b/(3*a**2*(a + b*x)**3) - b/(a**3*(a + b*x)**2) - 3*b/(a**4*(a + b*x)) - 1/(a**4*x) - 4*b*log(x)/a**5 + 4*b*log(a + b*x)/a**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_204():
    f = 1/(x**3*(a + b*x)**4)
    F = b**2/(3*a**3*(a + b*x)**3) + 3*b**2/(2*a**4*(a + b*x)**2) - 1/(2*a**4*x**2) + 6*b**2/(a**5*(a + b*x)) + 4*b/(a**5*x) + 10*b**2*log(x)/a**6 - 10*b**2*log(a + b*x)/a**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_205():
    f = 1/(x**4*(a + b*x)**4)
    F = -b**3/(3*a**4*(a + b*x)**3) - 1/(3*a**4*x**3) - 2*b**3/(a**5*(a + b*x)**2) + 2*b/(a**5*x**2) - 10*b**3/(a**6*(a + b*x)) - 10*b**2/(a**6*x) - 20*b**3*log(x)/a**7 + 20*b**3*log(a + b*x)/a**7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_206():
    f = 1/(x**5*(a + b*x)**4)
    F = -1/(4*a**4*x**4) + b**4/(3*a**5*(a + b*x)**3) + 4*b/(3*a**5*x**3) + 5*b**4/(2*a**6*(a + b*x)**2) - 5*b**2/(a**6*x**2) + 15*b**4/(a**7*(a + b*x)) + 20*b**3/(a**7*x) + 35*b**4*log(x)/a**8 - 35*b**4*log(a + b*x)/a**8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_207():
    f = x**10/(a + b*x)**7
    F = -a**10/(6*b**11*(a + b*x)**6) + 2*a**9/(b**11*(a + b*x)**5) - 45*a**8/(4*b**11*(a + b*x)**4) + 40*a**7/(b**11*(a + b*x)**3) - 105*a**6/(b**11*(a + b*x)**2) + 252*a**5/(b**11*(a + b*x)) + 210*a**4*log(a + b*x)/b**11 - 84*a**3*x/b**10 + 14*a**2*x**2/b**9 - 7*a*x**3/(3*b**8) + x**4/(4*b**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_208():
    f = x**9/(a + b*x)**7
    F = a**9/(6*b**10*(a + b*x)**6) - 9*a**8/(5*b**10*(a + b*x)**5) + 9*a**7/(b**10*(a + b*x)**4) - 28*a**6/(b**10*(a + b*x)**3) + 63*a**5/(b**10*(a + b*x)**2) - 126*a**4/(b**10*(a + b*x)) - 84*a**3*log(a + b*x)/b**10 + 28*a**2*x/b**9 - 7*a*x**2/(2*b**8) + x**3/(3*b**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_209():
    f = x**8/(a + b*x)**7
    F = -a**8/(6*b**9*(a + b*x)**6) + 8*a**7/(5*b**9*(a + b*x)**5) - 7*a**6/(b**9*(a + b*x)**4) + 56*a**5/(3*b**9*(a + b*x)**3) - 35*a**4/(b**9*(a + b*x)**2) + 56*a**3/(b**9*(a + b*x)) + 28*a**2*log(a + b*x)/b**9 - 7*a*x/b**8 + x**2/(2*b**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_210():
    f = x**7/(a + b*x)**7
    F = a**7/(6*b**8*(a + b*x)**6) - 7*a**6/(5*b**8*(a + b*x)**5) + 21*a**5/(4*b**8*(a + b*x)**4) - 35*a**4/(3*b**8*(a + b*x)**3) + 35*a**3/(2*b**8*(a + b*x)**2) - 21*a**2/(b**8*(a + b*x)) - 7*a*log(a + b*x)/b**8 + x/b**7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_211():
    f = x**6/(a + b*x)**7
    F = -a**6/(6*b**7*(a + b*x)**6) + 6*a**5/(5*b**7*(a + b*x)**5) - 15*a**4/(4*b**7*(a + b*x)**4) + 20*a**3/(3*b**7*(a + b*x)**3) - 15*a**2/(2*b**7*(a + b*x)**2) + 6*a/(b**7*(a + b*x)) + log(a + b*x)/b**7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_212():
    f = x**5/(a + b*x)**7
    F = x**6/(6*a*(a + b*x)**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_213():
    f = x**4/(a + b*x)**7
    F = x**5/(6*a*(a + b*x)**6) + x**5/(30*a**2*(a + b*x)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_214():
    f = x**2/(a + b*x)**7
    F = -a**2/(6*b**3*(a + b*x)**6) + 2*a/(5*b**3*(a + b*x)**5) - 1/(4*b**3*(a + b*x)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_215():
    f = x/(a + b*x)**7
    F = a/(6*b**2*(a + b*x)**6) - 1/(5*b**2*(a + b*x)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_216():
    f = (a + b*x)**(-7)
    F = -1/(6*b*(a + b*x)**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_217():
    f = 1/(x*(a + b*x)**7)
    F = 1/(6*a*(a + b*x)**6) + 1/(5*a**2*(a + b*x)**5) + 1/(4*a**3*(a + b*x)**4) + 1/(3*a**4*(a + b*x)**3) + 1/(2*a**5*(a + b*x)**2) + 1/(a**6*(a + b*x)) + log(x)/a**7 - log(a + b*x)/a**7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_218():
    f = 1/(x**2*(a + b*x)**7)
    F = -b/(6*a**2*(a + b*x)**6) - 2*b/(5*a**3*(a + b*x)**5) - 3*b/(4*a**4*(a + b*x)**4) - 4*b/(3*a**5*(a + b*x)**3) - 5*b/(2*a**6*(a + b*x)**2) - 6*b/(a**7*(a + b*x)) - 1/(a**7*x) - 7*b*log(x)/a**8 + 7*b*log(a + b*x)/a**8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_219():
    f = 1/(x**3*(a + b*x)**7)
    F = b**2/(6*a**3*(a + b*x)**6) + 3*b**2/(5*a**4*(a + b*x)**5) + 3*b**2/(2*a**5*(a + b*x)**4) + 10*b**2/(3*a**6*(a + b*x)**3) + 15*b**2/(2*a**7*(a + b*x)**2) - 1/(2*a**7*x**2) + 21*b**2/(a**8*(a + b*x)) + 7*b/(a**8*x) + 28*b**2*log(x)/a**9 - 28*b**2*log(a + b*x)/a**9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_220():
    f = 1/(x**4*(a + b*x)**7)
    F = -b**3/(6*a**4*(a + b*x)**6) - 4*b**3/(5*a**5*(a + b*x)**5) - 5*b**3/(2*a**6*(a + b*x)**4) - 20*b**3/(3*a**7*(a + b*x)**3) - 1/(3*a**7*x**3) - 35*b**3/(2*a**8*(a + b*x)**2) + 7*b/(2*a**8*x**2) - 56*b**3/(a**9*(a + b*x)) - 28*b**2/(a**9*x) - 84*b**3*log(x)/a**10 + 84*b**3*log(a + b*x)/a**10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_221():
    f = x**12/(a + b*x)**10
    F = -a**12/(9*b**13*(a + b*x)**9) + 3*a**11/(2*b**13*(a + b*x)**8) - 66*a**10/(7*b**13*(a + b*x)**7) + 110*a**9/(3*b**13*(a + b*x)**6) - 99*a**8/(b**13*(a + b*x)**5) + 198*a**7/(b**13*(a + b*x)**4) - 308*a**6/(b**13*(a + b*x)**3) + 396*a**5/(b**13*(a + b*x)**2) - 495*a**4/(b**13*(a + b*x)) - 220*a**3*log(a + b*x)/b**13 + 55*a**2*x/b**12 - 5*a*x**2/b**11 + x**3/(3*b**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_222():
    f = x**11/(a + b*x)**10
    F = a**11/(9*b**12*(a + b*x)**9) - 11*a**10/(8*b**12*(a + b*x)**8) + 55*a**9/(7*b**12*(a + b*x)**7) - 55*a**8/(2*b**12*(a + b*x)**6) + 66*a**7/(b**12*(a + b*x)**5) - 231*a**6/(2*b**12*(a + b*x)**4) + 154*a**5/(b**12*(a + b*x)**3) - 165*a**4/(b**12*(a + b*x)**2) + 165*a**3/(b**12*(a + b*x)) + 55*a**2*log(a + b*x)/b**12 - 10*a*x/b**11 + x**2/(2*b**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_223():
    f = x**10/(a + b*x)**10
    F = -a**10/(9*b**11*(a + b*x)**9) + 5*a**9/(4*b**11*(a + b*x)**8) - 45*a**8/(7*b**11*(a + b*x)**7) + 20*a**7/(b**11*(a + b*x)**6) - 42*a**6/(b**11*(a + b*x)**5) + 63*a**5/(b**11*(a + b*x)**4) - 70*a**4/(b**11*(a + b*x)**3) + 60*a**3/(b**11*(a + b*x)**2) - 45*a**2/(b**11*(a + b*x)) - 10*a*log(a + b*x)/b**11 + x/b**10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_224():
    f = x**9/(a + b*x)**10
    F = a**9/(9*b**10*(a + b*x)**9) - 9*a**8/(8*b**10*(a + b*x)**8) + 36*a**7/(7*b**10*(a + b*x)**7) - 14*a**6/(b**10*(a + b*x)**6) + 126*a**5/(5*b**10*(a + b*x)**5) - 63*a**4/(2*b**10*(a + b*x)**4) + 28*a**3/(b**10*(a + b*x)**3) - 18*a**2/(b**10*(a + b*x)**2) + 9*a/(b**10*(a + b*x)) + log(a + b*x)/b**10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_225():
    f = x**8/(a + b*x)**10
    F = x**9/(9*a*(a + b*x)**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_226():
    f = x**7/(a + b*x)**10
    F = x**8/(9*a*(a + b*x)**9) + x**8/(72*a**2*(a + b*x)**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_227():
    f = x**6/(a + b*x)**10
    F = x**7/(9*a*(a + b*x)**9) + x**7/(36*a**2*(a + b*x)**8) + x**7/(252*a**3*(a + b*x)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_228():
    f = x**5/(a + b*x)**10
    F = x**6/(9*a*(a + b*x)**9) + x**6/(24*a**2*(a + b*x)**8) + x**6/(84*a**3*(a + b*x)**7) + x**6/(504*a**4*(a + b*x)**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_229():
    f = x**4/(a + b*x)**10
    F = -a**4/(9*b**5*(a + b*x)**9) + a**3/(2*b**5*(a + b*x)**8) - 6*a**2/(7*b**5*(a + b*x)**7) + 2*a/(3*b**5*(a + b*x)**6) - 1/(5*b**5*(a + b*x)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_230():
    f = x**3/(a + b*x)**10
    F = a**3/(9*b**4*(a + b*x)**9) - 3*a**2/(8*b**4*(a + b*x)**8) + 3*a/(7*b**4*(a + b*x)**7) - 1/(6*b**4*(a + b*x)**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_231():
    f = x**2/(a + b*x)**10
    F = -a**2/(9*b**3*(a + b*x)**9) + a/(4*b**3*(a + b*x)**8) - 1/(7*b**3*(a + b*x)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_232():
    f = x/(a + b*x)**10
    F = a/(9*b**2*(a + b*x)**9) - 1/(8*b**2*(a + b*x)**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_233():
    f = (a + b*x)**(-10)
    F = -1/(9*b*(a + b*x)**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_234():
    f = 1/(x*(a + b*x)**10)
    F = 1/(9*a*(a + b*x)**9) + 1/(8*a**2*(a + b*x)**8) + 1/(7*a**3*(a + b*x)**7) + 1/(6*a**4*(a + b*x)**6) + 1/(5*a**5*(a + b*x)**5) + 1/(4*a**6*(a + b*x)**4) + 1/(3*a**7*(a + b*x)**3) + 1/(2*a**8*(a + b*x)**2) + 1/(a**9*(a + b*x)) + log(x)/a**10 - log(a + b*x)/a**10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_235():
    f = 1/(x**2*(a + b*x)**10)
    F = -b/(9*a**2*(a + b*x)**9) - b/(4*a**3*(a + b*x)**8) - 3*b/(7*a**4*(a + b*x)**7) - 2*b/(3*a**5*(a + b*x)**6) - b/(a**6*(a + b*x)**5) - 3*b/(2*a**7*(a + b*x)**4) - 7*b/(3*a**8*(a + b*x)**3) - 4*b/(a**9*(a + b*x)**2) - 9*b/(a**10*(a + b*x)) - 1/(a**10*x) - 10*b*log(x)/a**11 + 10*b*log(a + b*x)/a**11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_236():
    f = 1/(x**3*(a + b*x)**10)
    F = b**2/(9*a**3*(a + b*x)**9) + 3*b**2/(8*a**4*(a + b*x)**8) + 6*b**2/(7*a**5*(a + b*x)**7) + 5*b**2/(3*a**6*(a + b*x)**6) + 3*b**2/(a**7*(a + b*x)**5) + 21*b**2/(4*a**8*(a + b*x)**4) + 28*b**2/(3*a**9*(a + b*x)**3) + 18*b**2/(a**10*(a + b*x)**2) - 1/(2*a**10*x**2) + 45*b**2/(a**11*(a + b*x)) + 10*b/(a**11*x) + 55*b**2*log(x)/a**12 - 55*b**2*log(a + b*x)/a**12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_237():
    f = 1/(x**4*(a + b*x)**10)
    F = -b**3/(9*a**4*(a + b*x)**9) - b**3/(2*a**5*(a + b*x)**8) - 10*b**3/(7*a**6*(a + b*x)**7) - 10*b**3/(3*a**7*(a + b*x)**6) - 7*b**3/(a**8*(a + b*x)**5) - 14*b**3/(a**9*(a + b*x)**4) - 28*b**3/(a**10*(a + b*x)**3) - 1/(3*a**10*x**3) - 60*b**3/(a**11*(a + b*x)**2) + 5*b/(a**11*x**2) - 165*b**3/(a**12*(a + b*x)) - 55*b**2/(a**12*x) - 220*b**3*log(x)/a**13 + 220*b**3*log(a + b*x)/a**13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_238():
    f = (a + b*x)**12/x**10
    F = -a**12/(9*x**9) - 3*a**11*b/(2*x**8) - 66*a**10*b**2/(7*x**7) - 110*a**9*b**3/(3*x**6) - 99*a**8*b**4/x**5 - 198*a**7*b**5/x**4 - 308*a**6*b**6/x**3 - 396*a**5*b**7/x**2 - 495*a**4*b**8/x + 220*a**3*b**9*log(x) + 66*a**2*b**10*x + 6*a*b**11*x**2 + b**12*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_239():
    f = (a + b*x)**11/x**10
    F = -a**11/(9*x**9) - 11*a**10*b/(8*x**8) - 55*a**9*b**2/(7*x**7) - 55*a**8*b**3/(2*x**6) - 66*a**7*b**4/x**5 - 231*a**6*b**5/(2*x**4) - 154*a**5*b**6/x**3 - 165*a**4*b**7/x**2 - 165*a**3*b**8/x + 55*a**2*b**9*log(x) + 11*a*b**10*x + b**11*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_240():
    f = (a + b*x)**10/x**10
    F = -a**10/(9*x**9) - 5*a**9*b/(4*x**8) - 45*a**8*b**2/(7*x**7) - 20*a**7*b**3/x**6 - 42*a**6*b**4/x**5 - 63*a**5*b**5/x**4 - 70*a**4*b**6/x**3 - 60*a**3*b**7/x**2 - 45*a**2*b**8/x + 10*a*b**9*log(x) + b**10*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_241():
    f = (a + b*x)**9/x**10
    F = -a**9/(9*x**9) - 9*a**8*b/(8*x**8) - 36*a**7*b**2/(7*x**7) - 14*a**6*b**3/x**6 - 126*a**5*b**4/(5*x**5) - 63*a**4*b**5/(2*x**4) - 28*a**3*b**6/x**3 - 18*a**2*b**7/x**2 - 9*a*b**8/x + b**9*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_242():
    f = (a + b*x)**8/x**10
    F = -(a + b*x)**9/(9*a*x**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_243():
    f = (a + b*x)**7/x**10
    F = -(a + b*x)**8/(9*a*x**9) + b*(a + b*x)**8/(72*a**2*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_244():
    f = (a + b*x)**6/x**10
    F = -(a + b*x)**7/(9*a*x**9) + b*(a + b*x)**7/(36*a**2*x**8) - b**2*(a + b*x)**7/(252*a**3*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_245():
    f = (a + b*x)**5/x**10
    F = -a**5/(9*x**9) - 5*a**4*b/(8*x**8) - 10*a**3*b**2/(7*x**7) - 5*a**2*b**3/(3*x**6) - a*b**4/x**5 - b**5/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_246():
    f = (a + b*x)**4/x**10
    F = -a**4/(9*x**9) - a**3*b/(2*x**8) - 6*a**2*b**2/(7*x**7) - 2*a*b**3/(3*x**6) - b**4/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_247():
    f = (a + b*x)**3/x**10
    F = -a**3/(9*x**9) - 3*a**2*b/(8*x**8) - 3*a*b**2/(7*x**7) - b**3/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_248():
    f = (a + b*x)**2/x**10
    F = -a**2/(9*x**9) - a*b/(4*x**8) - b**2/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_249():
    f = (a + b*x)/x**10
    F = -a/(9*x**9) - b/(8*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_250():
    f = x**(-10)
    F = -1/(9*x**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_251():
    f = 1/(x**10*(a + b*x))
    F = -1/(9*a*x**9) + b/(8*a**2*x**8) - b**2/(7*a**3*x**7) + b**3/(6*a**4*x**6) - b**4/(5*a**5*x**5) + b**5/(4*a**6*x**4) - b**6/(3*a**7*x**3) + b**7/(2*a**8*x**2) - b**8/(a**9*x) - b**9*log(x)/a**10 + b**9*log(a + b*x)/a**10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_252():
    f = 1/(x**10*(a + b*x)**2)
    F = -1/(9*a**2*x**9) + b/(4*a**3*x**8) - 3*b**2/(7*a**4*x**7) + 2*b**3/(3*a**5*x**6) - b**4/(a**6*x**5) + 3*b**5/(2*a**7*x**4) - 7*b**6/(3*a**8*x**3) + 4*b**7/(a**9*x**2) - b**9/(a**10*(a + b*x)) - 9*b**8/(a**10*x) - 10*b**9*log(x)/a**11 + 10*b**9*log(a + b*x)/a**11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_253():
    f = 1/(x**10*(a + b*x)**3)
    F = -1/(9*a**3*x**9) + 3*b/(8*a**4*x**8) - 6*b**2/(7*a**5*x**7) + 5*b**3/(3*a**6*x**6) - 3*b**4/(a**7*x**5) + 21*b**5/(4*a**8*x**4) - 28*b**6/(3*a**9*x**3) - b**9/(2*a**10*(a + b*x)**2) + 18*b**7/(a**10*x**2) - 10*b**9/(a**11*(a + b*x)) - 45*b**8/(a**11*x) - 55*b**9*log(x)/a**12 + 55*b**9*log(a + b*x)/a**12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_254():
    f = 1/(x*(3*x + 2))
    F = log(x)/2 - log(3*x + 2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_255():
    f = 1/(x*(6*x + 4))
    F = log(x)/4 - log(3*x + 2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_256():
    f = 1/(x**2*(6*x + 4))
    F = -3*log(x)/8 + 3*log(3*x + 2)/8 - 1/(4*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_257():
    f = 1/(x**3*(6*x + 4))
    F = 9*log(x)/16 - 9*log(3*x + 2)/16 + 3/(8*x) - 1/(8*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_258():
    f = 1/(x**4*(6*x + 4))
    F = -27*log(x)/32 + 27*log(3*x + 2)/32 - 9/(16*x) + 3/(16*x**2) - 1/(12*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_259():
    f = 1/(x**5*(6*x + 4))
    F = 81*log(x)/64 - 81*log(3*x + 2)/64 + 27/(32*x) - 9/(32*x**2) + 1/(8*x**3) - 1/(16*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_260():
    f = 1/(x*(6*x + 4)**2)
    F = log(x)/16 - log(3*x + 2)/16 + 1/(24*x + 16)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_261():
    f = 1/(x**2*(6*x + 4)**2)
    F = -3*log(x)/16 + 3*log(3*x + 2)/16 - 3/(48*x + 32) - 1/(16*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_262():
    f = 1/(x**3*(6*x + 4)**2)
    F = 27*log(x)/64 - 27*log(3*x + 2)/64 + 9/(96*x + 64) + 3/(16*x) - 1/(32*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_263():
    f = 1/(x**4*(6*x + 4)**2)
    F = -27*log(x)/32 + 27*log(3*x + 2)/32 - 27/(192*x + 128) - 27/(64*x) + 3/(32*x**2) - 1/(48*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_264():
    f = 1/(x**5*(6*x + 4)**2)
    F = 405*log(x)/256 - 405*log(3*x + 2)/256 + 81/(384*x + 256) + 27/(32*x) - 27/(128*x**2) + 1/(16*x**3) - 1/(64*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_265():
    f = 1/(x*(6*x + 4)**3)
    F = log(x)/64 - log(3*x + 2)/64 + 1/(96*x + 64) + 1/(32*(3*x + 2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_266():
    f = 1/(x**2*(6*x + 4)**3)
    F = -9*log(x)/128 + 9*log(3*x + 2)/128 - 3/(96*x + 64) - 3/(64*(3*x + 2)**2) - 1/(64*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_267():
    f = 1/(x**3*(6*x + 4)**3)
    F = 27*log(x)/128 - 27*log(3*x + 2)/128 + 27/(384*x + 256) + 9/(128*(3*x + 2)**2) + 9/(128*x) - 1/(128*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_268():
    f = 1/(x**4*(6*x + 4)**3)
    F = -135*log(x)/256 + 135*log(3*x + 2)/256 - 27/(192*x + 128) - 27/(256*(3*x + 2)**2) - 27/(128*x) + 9/(256*x**2) - 1/(192*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_269():
    f = 1/(x**5*(6*x + 4)**3)
    F = 1215*log(x)/1024 - 1215*log(3*x + 2)/1024 + 405/(1536*x + 1024) + 81/(512*(3*x + 2)**2) + 135/(256*x) - 27/(256*x**2) + 3/(128*x**3) - 1/(256*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_270():
    f = 1/(2*x + 2)
    F = log(x + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_271():
    f = 1/(4 - 6*x)
    F = -log(2 - 3*x)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_272():
    f = 1/(sqrt(a)*x + a)
    F = log(sqrt(a) + x)/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_273():
    f = 1/(a + x*sqrt(-a))
    F = log(a + x*sqrt(-a))/sqrt(-a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_274():
    f = 1/(a**2 + x*sqrt(-a))
    F = log(a**2 + x*sqrt(-a))/sqrt(-a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_275():
    f = 1/(a**3 + x*sqrt(-a))
    F = log(a**3 + x*sqrt(-a))/sqrt(-a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_276():
    f = 1/(x*sqrt(-a) + 1/a)
    F = log(-x*(-a)**(sympy.S(3)/2) + 1)/sqrt(-a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_277():
    f = 1/(x*sqrt(-a) + a**(-2))
    F = log(x*(-a)**(sympy.S(5)/2) + 1)/sqrt(-a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_278():
    f = 1/(x*(b*x + 1))
    F = log(x) - log(b*x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_279():
    f = 1/(x*(b*x - 1))
    F = -log(x) + log(-b*x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_280():
    f = 1/(x**2*(b*x + 1))
    F = -b*log(x) + b*log(b*x + 1) - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_281():
    f = 1/(x**2*(b*x - 1))
    F = -b*log(x) + b*log(-b*x + 1) + 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_282():
    f = b/x + 1/(x**2*(b*x + 1))
    F = b*log(b*x + 1) - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_283():
    f = x**3*sqrt(a + b*x)
    F = -2*a**3*(a + b*x)**(sympy.S(3)/2)/(3*b**4) + 6*a**2*(a + b*x)**(sympy.S(5)/2)/(5*b**4) - 6*a*(a + b*x)**(sympy.S(7)/2)/(7*b**4) + 2*(a + b*x)**(sympy.S(9)/2)/(9*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_284():
    f = x**2*sqrt(a + b*x)
    F = 2*a**2*(a + b*x)**(sympy.S(3)/2)/(3*b**3) - 4*a*(a + b*x)**(sympy.S(5)/2)/(5*b**3) + 2*(a + b*x)**(sympy.S(7)/2)/(7*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_285():
    f = x*sqrt(a + b*x)
    F = -2*a*(a + b*x)**(sympy.S(3)/2)/(3*b**2) + 2*(a + b*x)**(sympy.S(5)/2)/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_286():
    f = sqrt(a + b*x)
    F = 2*(a + b*x)**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_287():
    f = sqrt(a + b*x)/x
    F = -2*sqrt(a)*atanh(sqrt(a + b*x)/sqrt(a)) + 2*sqrt(a + b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_288():
    f = sqrt(a + b*x)/x**2
    F = -sqrt(a + b*x)/x - b*atanh(sqrt(a + b*x)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_289():
    f = sqrt(a + b*x)/x**3
    F = -sqrt(a + b*x)/(2*x**2) - b*sqrt(a + b*x)/(4*a*x) + b**2*atanh(sqrt(a + b*x)/sqrt(a))/(4*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_290():
    f = sqrt(a + b*x)/x**4
    F = -sqrt(a + b*x)/(3*x**3) - b*sqrt(a + b*x)/(12*a*x**2) + b**2*sqrt(a + b*x)/(8*a**2*x) - b**3*atanh(sqrt(a + b*x)/sqrt(a))/(8*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_291():
    f = x**3*(a + b*x)**(sympy.S(3)/2)
    F = -2*a**3*(a + b*x)**(sympy.S(5)/2)/(5*b**4) + 6*a**2*(a + b*x)**(sympy.S(7)/2)/(7*b**4) - 2*a*(a + b*x)**(sympy.S(9)/2)/(3*b**4) + 2*(a + b*x)**(sympy.S(11)/2)/(11*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_292():
    f = x**2*(a + b*x)**(sympy.S(3)/2)
    F = 2*a**2*(a + b*x)**(sympy.S(5)/2)/(5*b**3) - 4*a*(a + b*x)**(sympy.S(7)/2)/(7*b**3) + 2*(a + b*x)**(sympy.S(9)/2)/(9*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_293():
    f = x*(a + b*x)**(sympy.S(3)/2)
    F = -2*a*(a + b*x)**(sympy.S(5)/2)/(5*b**2) + 2*(a + b*x)**(sympy.S(7)/2)/(7*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_294():
    f = (a + b*x)**(sympy.S(3)/2)
    F = 2*(a + b*x)**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_295():
    f = (a + b*x)**(sympy.S(3)/2)/x
    F = -2*a**(sympy.S(3)/2)*atanh(sqrt(a + b*x)/sqrt(a)) + 2*a*sqrt(a + b*x) + 2*(a + b*x)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_296():
    f = (a + b*x)**(sympy.S(3)/2)/x**2
    F = -3*sqrt(a)*b*atanh(sqrt(a + b*x)/sqrt(a)) + 3*b*sqrt(a + b*x) - (a + b*x)**(sympy.S(3)/2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_297():
    f = (a + b*x)**(sympy.S(3)/2)/x**3
    F = -3*b*sqrt(a + b*x)/(4*x) - (a + b*x)**(sympy.S(3)/2)/(2*x**2) - 3*b**2*atanh(sqrt(a + b*x)/sqrt(a))/(4*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_298():
    f = (a + b*x)**(sympy.S(3)/2)/x**4
    F = -b*sqrt(a + b*x)/(4*x**2) - (a + b*x)**(sympy.S(3)/2)/(3*x**3) - b**2*sqrt(a + b*x)/(8*a*x) + b**3*atanh(sqrt(a + b*x)/sqrt(a))/(8*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_299():
    f = x**3*(a + b*x)**(sympy.S(5)/2)
    F = -2*a**3*(a + b*x)**(sympy.S(7)/2)/(7*b**4) + 2*a**2*(a + b*x)**(sympy.S(9)/2)/(3*b**4) - 6*a*(a + b*x)**(sympy.S(11)/2)/(11*b**4) + 2*(a + b*x)**(sympy.S(13)/2)/(13*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_300():
    f = x**2*(a + b*x)**(sympy.S(5)/2)
    F = 2*a**2*(a + b*x)**(sympy.S(7)/2)/(7*b**3) - 4*a*(a + b*x)**(sympy.S(9)/2)/(9*b**3) + 2*(a + b*x)**(sympy.S(11)/2)/(11*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_301():
    f = x*(a + b*x)**(sympy.S(5)/2)
    F = -2*a*(a + b*x)**(sympy.S(7)/2)/(7*b**2) + 2*(a + b*x)**(sympy.S(9)/2)/(9*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_302():
    f = (a + b*x)**(sympy.S(5)/2)
    F = 2*(a + b*x)**(sympy.S(7)/2)/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_303():
    f = (a + b*x)**(sympy.S(5)/2)/x
    F = -2*a**(sympy.S(5)/2)*atanh(sqrt(a + b*x)/sqrt(a)) + 2*a**2*sqrt(a + b*x) + 2*a*(a + b*x)**(sympy.S(3)/2)/3 + 2*(a + b*x)**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_304():
    f = (a + b*x)**(sympy.S(5)/2)/x**2
    F = -5*a**(sympy.S(3)/2)*b*atanh(sqrt(a + b*x)/sqrt(a)) + 5*a*b*sqrt(a + b*x) + 5*b*(a + b*x)**(sympy.S(3)/2)/3 - (a + b*x)**(sympy.S(5)/2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_305():
    f = (a + b*x)**(sympy.S(5)/2)/x**3
    F = -15*sqrt(a)*b**2*atanh(sqrt(a + b*x)/sqrt(a))/4 + 15*b**2*sqrt(a + b*x)/4 - 5*b*(a + b*x)**(sympy.S(3)/2)/(4*x) - (a + b*x)**(sympy.S(5)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_306():
    f = (a + b*x)**(sympy.S(5)/2)/x**4
    F = -5*b**2*sqrt(a + b*x)/(8*x) - 5*b*(a + b*x)**(sympy.S(3)/2)/(12*x**2) - (a + b*x)**(sympy.S(5)/2)/(3*x**3) - 5*b**3*atanh(sqrt(a + b*x)/sqrt(a))/(8*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_307():
    f = (a + b*x)**(sympy.S(5)/2)/x**5
    F = -5*b**2*sqrt(a + b*x)/(32*x**2) - 5*b*(a + b*x)**(sympy.S(3)/2)/(24*x**3) - (a + b*x)**(sympy.S(5)/2)/(4*x**4) - 5*b**3*sqrt(a + b*x)/(64*a*x) + 5*b**4*atanh(sqrt(a + b*x)/sqrt(a))/(64*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_308():
    f = x**7*(a + b*x)**(sympy.S(9)/2)
    F = -2*a**7*(a + b*x)**(sympy.S(11)/2)/(11*b**8) + 14*a**6*(a + b*x)**(sympy.S(13)/2)/(13*b**8) - 14*a**5*(a + b*x)**(sympy.S(15)/2)/(5*b**8) + 70*a**4*(a + b*x)**(sympy.S(17)/2)/(17*b**8) - 70*a**3*(a + b*x)**(sympy.S(19)/2)/(19*b**8) + 2*a**2*(a + b*x)**(sympy.S(21)/2)/b**8 - 14*a*(a + b*x)**(sympy.S(23)/2)/(23*b**8) + 2*(a + b*x)**(sympy.S(25)/2)/(25*b**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_309():
    f = x**6*(a + b*x)**(sympy.S(9)/2)
    F = 2*a**6*(a + b*x)**(sympy.S(11)/2)/(11*b**7) - 12*a**5*(a + b*x)**(sympy.S(13)/2)/(13*b**7) + 2*a**4*(a + b*x)**(sympy.S(15)/2)/b**7 - 40*a**3*(a + b*x)**(sympy.S(17)/2)/(17*b**7) + 30*a**2*(a + b*x)**(sympy.S(19)/2)/(19*b**7) - 4*a*(a + b*x)**(sympy.S(21)/2)/(7*b**7) + 2*(a + b*x)**(sympy.S(23)/2)/(23*b**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_310():
    f = x**5*(a + b*x)**(sympy.S(9)/2)
    F = -2*a**5*(a + b*x)**(sympy.S(11)/2)/(11*b**6) + 10*a**4*(a + b*x)**(sympy.S(13)/2)/(13*b**6) - 4*a**3*(a + b*x)**(sympy.S(15)/2)/(3*b**6) + 20*a**2*(a + b*x)**(sympy.S(17)/2)/(17*b**6) - 10*a*(a + b*x)**(sympy.S(19)/2)/(19*b**6) + 2*(a + b*x)**(sympy.S(21)/2)/(21*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_311():
    f = x**4*(a + b*x)**(sympy.S(9)/2)
    F = 2*a**4*(a + b*x)**(sympy.S(11)/2)/(11*b**5) - 8*a**3*(a + b*x)**(sympy.S(13)/2)/(13*b**5) + 4*a**2*(a + b*x)**(sympy.S(15)/2)/(5*b**5) - 8*a*(a + b*x)**(sympy.S(17)/2)/(17*b**5) + 2*(a + b*x)**(sympy.S(19)/2)/(19*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_312():
    f = x**3*(a + b*x)**(sympy.S(9)/2)
    F = -2*a**3*(a + b*x)**(sympy.S(11)/2)/(11*b**4) + 6*a**2*(a + b*x)**(sympy.S(13)/2)/(13*b**4) - 2*a*(a + b*x)**(sympy.S(15)/2)/(5*b**4) + 2*(a + b*x)**(sympy.S(17)/2)/(17*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_313():
    f = x**2*(a + b*x)**(sympy.S(9)/2)
    F = 2*a**2*(a + b*x)**(sympy.S(11)/2)/(11*b**3) - 4*a*(a + b*x)**(sympy.S(13)/2)/(13*b**3) + 2*(a + b*x)**(sympy.S(15)/2)/(15*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_314():
    f = x*(a + b*x)**(sympy.S(9)/2)
    F = -2*a*(a + b*x)**(sympy.S(11)/2)/(11*b**2) + 2*(a + b*x)**(sympy.S(13)/2)/(13*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_315():
    f = (a + b*x)**(sympy.S(9)/2)
    F = 2*(a + b*x)**(sympy.S(11)/2)/(11*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_316():
    f = (a + b*x)**(sympy.S(9)/2)/x
    F = -2*a**(sympy.S(9)/2)*atanh(sqrt(a + b*x)/sqrt(a)) + 2*a**4*sqrt(a + b*x) + 2*a**3*(a + b*x)**(sympy.S(3)/2)/3 + 2*a**2*(a + b*x)**(sympy.S(5)/2)/5 + 2*a*(a + b*x)**(sympy.S(7)/2)/7 + 2*(a + b*x)**(sympy.S(9)/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_317():
    f = (a + b*x)**(sympy.S(9)/2)/x**2
    F = -9*a**(sympy.S(7)/2)*b*atanh(sqrt(a + b*x)/sqrt(a)) + 9*a**3*b*sqrt(a + b*x) + 3*a**2*b*(a + b*x)**(sympy.S(3)/2) + 9*a*b*(a + b*x)**(sympy.S(5)/2)/5 + 9*b*(a + b*x)**(sympy.S(7)/2)/7 - (a + b*x)**(sympy.S(9)/2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_318():
    f = (a + b*x)**(sympy.S(9)/2)/x**3
    F = -63*a**(sympy.S(5)/2)*b**2*atanh(sqrt(a + b*x)/sqrt(a))/4 + 63*a**2*b**2*sqrt(a + b*x)/4 + 21*a*b**2*(a + b*x)**(sympy.S(3)/2)/4 + 63*b**2*(a + b*x)**(sympy.S(5)/2)/20 - 9*b*(a + b*x)**(sympy.S(7)/2)/(4*x) - (a + b*x)**(sympy.S(9)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_319():
    f = (a + b*x)**(sympy.S(9)/2)/x**4
    F = -105*a**(sympy.S(3)/2)*b**3*atanh(sqrt(a + b*x)/sqrt(a))/8 + 105*a*b**3*sqrt(a + b*x)/8 + 35*b**3*(a + b*x)**(sympy.S(3)/2)/8 - 21*b**2*(a + b*x)**(sympy.S(5)/2)/(8*x) - 3*b*(a + b*x)**(sympy.S(7)/2)/(4*x**2) - (a + b*x)**(sympy.S(9)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_320():
    f = (a + b*x)**(sympy.S(9)/2)/x**5
    F = -315*sqrt(a)*b**4*atanh(sqrt(a + b*x)/sqrt(a))/64 + 315*b**4*sqrt(a + b*x)/64 - 105*b**3*(a + b*x)**(sympy.S(3)/2)/(64*x) - 21*b**2*(a + b*x)**(sympy.S(5)/2)/(32*x**2) - 3*b*(a + b*x)**(sympy.S(7)/2)/(8*x**3) - (a + b*x)**(sympy.S(9)/2)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_321():
    f = (a + b*x)**(sympy.S(9)/2)/x**6
    F = -63*b**4*sqrt(a + b*x)/(128*x) - 21*b**3*(a + b*x)**(sympy.S(3)/2)/(64*x**2) - 21*b**2*(a + b*x)**(sympy.S(5)/2)/(80*x**3) - 9*b*(a + b*x)**(sympy.S(7)/2)/(40*x**4) - (a + b*x)**(sympy.S(9)/2)/(5*x**5) - 63*b**5*atanh(sqrt(a + b*x)/sqrt(a))/(128*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_322():
    f = (a + b*x)**(sympy.S(9)/2)/x**7
    F = -21*b**4*sqrt(a + b*x)/(256*x**2) - 7*b**3*(a + b*x)**(sympy.S(3)/2)/(64*x**3) - 21*b**2*(a + b*x)**(sympy.S(5)/2)/(160*x**4) - 3*b*(a + b*x)**(sympy.S(7)/2)/(20*x**5) - (a + b*x)**(sympy.S(9)/2)/(6*x**6) - 21*b**5*sqrt(a + b*x)/(512*a*x) + 21*b**6*atanh(sqrt(a + b*x)/sqrt(a))/(512*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_323():
    f = (a + b*x)**(sympy.S(9)/2)/x**8
    F = -3*b**4*sqrt(a + b*x)/(128*x**3) - 3*b**3*(a + b*x)**(sympy.S(3)/2)/(64*x**4) - 3*b**2*(a + b*x)**(sympy.S(5)/2)/(40*x**5) - 3*b*(a + b*x)**(sympy.S(7)/2)/(28*x**6) - (a + b*x)**(sympy.S(9)/2)/(7*x**7) - 3*b**5*sqrt(a + b*x)/(512*a*x**2) + 9*b**6*sqrt(a + b*x)/(1024*a**2*x) - 9*b**7*atanh(sqrt(a + b*x)/sqrt(a))/(1024*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_324():
    f = sqrt(-a + b*x)/x
    F = -2*sqrt(a)*atan(sqrt(-a + b*x)/sqrt(a)) + 2*sqrt(-a + b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_325():
    f = sqrt(-a + b*x)/x**2
    F = -sqrt(-a + b*x)/x + b*atan(sqrt(-a + b*x)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_326():
    f = sqrt(-a + b*x)/x**3
    F = -sqrt(-a + b*x)/(2*x**2) + b*sqrt(-a + b*x)/(4*a*x) + b**2*atan(sqrt(-a + b*x)/sqrt(a))/(4*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_327():
    f = (-a + b*x)**(sympy.S(3)/2)/x
    F = 2*a**(sympy.S(3)/2)*atan(sqrt(-a + b*x)/sqrt(a)) - 2*a*sqrt(-a + b*x) + 2*(-a + b*x)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_328():
    f = (-a + b*x)**(sympy.S(3)/2)/x**2
    F = -3*sqrt(a)*b*atan(sqrt(-a + b*x)/sqrt(a)) + 3*b*sqrt(-a + b*x) - (-a + b*x)**(sympy.S(3)/2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_329():
    f = (-a + b*x)**(sympy.S(3)/2)/x**3
    F = -3*b*sqrt(-a + b*x)/(4*x) - (-a + b*x)**(sympy.S(3)/2)/(2*x**2) + 3*b**2*atan(sqrt(-a + b*x)/sqrt(a))/(4*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_330():
    f = (-a + b*x)**(sympy.S(5)/2)/x
    F = -2*a**(sympy.S(5)/2)*atan(sqrt(-a + b*x)/sqrt(a)) + 2*a**2*sqrt(-a + b*x) - 2*a*(-a + b*x)**(sympy.S(3)/2)/3 + 2*(-a + b*x)**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_331():
    f = (-a + b*x)**(sympy.S(5)/2)/x**2
    F = 5*a**(sympy.S(3)/2)*b*atan(sqrt(-a + b*x)/sqrt(a)) - 5*a*b*sqrt(-a + b*x) + 5*b*(-a + b*x)**(sympy.S(3)/2)/3 - (-a + b*x)**(sympy.S(5)/2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_332():
    f = (-a + b*x)**(sympy.S(5)/2)/x**3
    F = -15*sqrt(a)*b**2*atan(sqrt(-a + b*x)/sqrt(a))/4 + 15*b**2*sqrt(-a + b*x)/4 - 5*b*(-a + b*x)**(sympy.S(3)/2)/(4*x) - (-a + b*x)**(sympy.S(5)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_333():
    f = x**4/sqrt(a + b*x)
    F = 2*a**4*sqrt(a + b*x)/b**5 - 8*a**3*(a + b*x)**(sympy.S(3)/2)/(3*b**5) + 12*a**2*(a + b*x)**(sympy.S(5)/2)/(5*b**5) - 8*a*(a + b*x)**(sympy.S(7)/2)/(7*b**5) + 2*(a + b*x)**(sympy.S(9)/2)/(9*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_334():
    f = x**3/sqrt(a + b*x)
    F = -2*a**3*sqrt(a + b*x)/b**4 + 2*a**2*(a + b*x)**(sympy.S(3)/2)/b**4 - 6*a*(a + b*x)**(sympy.S(5)/2)/(5*b**4) + 2*(a + b*x)**(sympy.S(7)/2)/(7*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_335():
    f = x**2/sqrt(a + b*x)
    F = 2*a**2*sqrt(a + b*x)/b**3 - 4*a*(a + b*x)**(sympy.S(3)/2)/(3*b**3) + 2*(a + b*x)**(sympy.S(5)/2)/(5*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_336():
    f = x/sqrt(a + b*x)
    F = -2*a*sqrt(a + b*x)/b**2 + 2*(a + b*x)**(sympy.S(3)/2)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_337():
    f = 1/sqrt(a + b*x)
    F = 2*sqrt(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_338():
    f = 1/(x*sqrt(a + b*x))
    F = -2*atanh(sqrt(a + b*x)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_339():
    f = 1/(x**2*sqrt(a + b*x))
    F = -sqrt(a + b*x)/(a*x) + b*atanh(sqrt(a + b*x)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_340():
    f = 1/(x**3*sqrt(a + b*x))
    F = -sqrt(a + b*x)/(2*a*x**2) + 3*b*sqrt(a + b*x)/(4*a**2*x) - 3*b**2*atanh(sqrt(a + b*x)/sqrt(a))/(4*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_341():
    f = 1/(x**4*sqrt(a + b*x))
    F = -sqrt(a + b*x)/(3*a*x**3) + 5*b*sqrt(a + b*x)/(12*a**2*x**2) - 5*b**2*sqrt(a + b*x)/(8*a**3*x) + 5*b**3*atanh(sqrt(a + b*x)/sqrt(a))/(8*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_342():
    f = x**4/(a + b*x)**(sympy.S(3)/2)
    F = -2*a**4/(b**5*sqrt(a + b*x)) - 8*a**3*sqrt(a + b*x)/b**5 + 4*a**2*(a + b*x)**(sympy.S(3)/2)/b**5 - 8*a*(a + b*x)**(sympy.S(5)/2)/(5*b**5) + 2*(a + b*x)**(sympy.S(7)/2)/(7*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_343():
    f = x**3/(a + b*x)**(sympy.S(3)/2)
    F = 2*a**3/(b**4*sqrt(a + b*x)) + 6*a**2*sqrt(a + b*x)/b**4 - 2*a*(a + b*x)**(sympy.S(3)/2)/b**4 + 2*(a + b*x)**(sympy.S(5)/2)/(5*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_344():
    f = x**2/(a + b*x)**(sympy.S(3)/2)
    F = -2*a**2/(b**3*sqrt(a + b*x)) - 4*a*sqrt(a + b*x)/b**3 + 2*(a + b*x)**(sympy.S(3)/2)/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_345():
    f = x/(a + b*x)**(sympy.S(3)/2)
    F = 2*a/(b**2*sqrt(a + b*x)) + 2*sqrt(a + b*x)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_346():
    f = (a + b*x)**(sympy.S(-3)/2)
    F = -2/(b*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_347():
    f = 1/(x*(a + b*x)**(sympy.S(3)/2))
    F = 2/(a*sqrt(a + b*x)) - 2*atanh(sqrt(a + b*x)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_348():
    f = 1/(x**2*(a + b*x)**(sympy.S(3)/2))
    F = -1/(a*x*sqrt(a + b*x)) - 3*b/(a**2*sqrt(a + b*x)) + 3*b*atanh(sqrt(a + b*x)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_349():
    f = 1/(x**3*(a + b*x)**(sympy.S(3)/2))
    F = -1/(2*a*x**2*sqrt(a + b*x)) + 5*b/(4*a**2*x*sqrt(a + b*x)) + 15*b**2/(4*a**3*sqrt(a + b*x)) - 15*b**2*atanh(sqrt(a + b*x)/sqrt(a))/(4*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_350():
    f = x**4/(a + b*x)**(sympy.S(5)/2)
    F = -2*a**4/(3*b**5*(a + b*x)**(sympy.S(3)/2)) + 8*a**3/(b**5*sqrt(a + b*x)) + 12*a**2*sqrt(a + b*x)/b**5 - 8*a*(a + b*x)**(sympy.S(3)/2)/(3*b**5) + 2*(a + b*x)**(sympy.S(5)/2)/(5*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_351():
    f = x**3/(a + b*x)**(sympy.S(5)/2)
    F = 2*a**3/(3*b**4*(a + b*x)**(sympy.S(3)/2)) - 6*a**2/(b**4*sqrt(a + b*x)) - 6*a*sqrt(a + b*x)/b**4 + 2*(a + b*x)**(sympy.S(3)/2)/(3*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_352():
    f = x**2/(a + b*x)**(sympy.S(5)/2)
    F = -2*a**2/(3*b**3*(a + b*x)**(sympy.S(3)/2)) + 4*a/(b**3*sqrt(a + b*x)) + 2*sqrt(a + b*x)/b**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_353():
    f = x/(a + b*x)**(sympy.S(5)/2)
    F = 2*a/(3*b**2*(a + b*x)**(sympy.S(3)/2)) - 2/(b**2*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_354():
    f = (a + b*x)**(sympy.S(-5)/2)
    F = -2/(3*b*(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_355():
    f = 1/(x*(a + b*x)**(sympy.S(5)/2))
    F = 2/(3*a*(a + b*x)**(sympy.S(3)/2)) + 2/(a**2*sqrt(a + b*x)) - 2*atanh(sqrt(a + b*x)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_356():
    f = 1/(x**2*(a + b*x)**(sympy.S(5)/2))
    F = -1/(a*x*(a + b*x)**(sympy.S(3)/2)) - 5*b/(3*a**2*(a + b*x)**(sympy.S(3)/2)) - 5*b/(a**3*sqrt(a + b*x)) + 5*b*atanh(sqrt(a + b*x)/sqrt(a))/a**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_357():
    f = 1/(x**3*(a + b*x)**(sympy.S(5)/2))
    F = -1/(2*a*x**2*(a + b*x)**(sympy.S(3)/2)) + 7*b/(4*a**2*x*(a + b*x)**(sympy.S(3)/2)) + 35*b**2/(12*a**3*(a + b*x)**(sympy.S(3)/2)) + 35*b**2/(4*a**4*sqrt(a + b*x)) - 35*b**2*atanh(sqrt(a + b*x)/sqrt(a))/(4*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_358():
    f = 1/(x*sqrt(-a + b*x))
    F = 2*atan(sqrt(-a + b*x)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_359():
    f = 1/(x**2*sqrt(-a + b*x))
    F = sqrt(-a + b*x)/(a*x) + b*atan(sqrt(-a + b*x)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_360():
    f = 1/(x**3*sqrt(-a + b*x))
    F = sqrt(-a + b*x)/(2*a*x**2) + 3*b*sqrt(-a + b*x)/(4*a**2*x) + 3*b**2*atan(sqrt(-a + b*x)/sqrt(a))/(4*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_361():
    f = 1/(x*(-a + b*x)**(sympy.S(3)/2))
    F = -2/(a*sqrt(-a + b*x)) - 2*atan(sqrt(-a + b*x)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_362():
    f = 1/(x**2*(-a + b*x)**(sympy.S(3)/2))
    F = 1/(a*x*sqrt(-a + b*x)) - 3*b/(a**2*sqrt(-a + b*x)) - 3*b*atan(sqrt(-a + b*x)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_363():
    f = 1/(x**3*(-a + b*x)**(sympy.S(3)/2))
    F = 1/(2*a*x**2*sqrt(-a + b*x)) + 5*b/(4*a**2*x*sqrt(-a + b*x)) - 15*b**2/(4*a**3*sqrt(-a + b*x)) - 15*b**2*atan(sqrt(-a + b*x)/sqrt(a))/(4*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_364():
    f = 1/(x*(-a + b*x)**(sympy.S(5)/2))
    F = -2/(3*a*(-a + b*x)**(sympy.S(3)/2)) + 2/(a**2*sqrt(-a + b*x)) + 2*atan(sqrt(-a + b*x)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_365():
    f = 1/(x**2*(-a + b*x)**(sympy.S(5)/2))
    F = 1/(a*x*(-a + b*x)**(sympy.S(3)/2)) - 5*b/(3*a**2*(-a + b*x)**(sympy.S(3)/2)) + 5*b/(a**3*sqrt(-a + b*x)) + 5*b*atan(sqrt(-a + b*x)/sqrt(a))/a**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_366():
    f = 1/(x**3*(-a + b*x)**(sympy.S(5)/2))
    F = 1/(2*a*x**2*(-a + b*x)**(sympy.S(3)/2)) + 7*b/(4*a**2*x*(-a + b*x)**(sympy.S(3)/2)) - 35*b**2/(12*a**3*(-a + b*x)**(sympy.S(3)/2)) + 35*b**2/(4*a**4*sqrt(-a + b*x)) + 35*b**2*atan(sqrt(-a + b*x)/sqrt(a))/(4*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_367():
    f = x**(m - 1)*(2*a*m + b*x*(2*m - 1))/(2*(a + b*x)**(sympy.S(3)/2))
    F = x**m/sqrt(a + b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_368():
    f = -b*x**m/(2*(a + b*x)**(sympy.S(3)/2)) + m*x**(m - 1)/sqrt(a + b*x)
    F = x**m/sqrt(a + b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_369():
    f = 1/(x*sqrt(a + b*x))
    F = -2*atanh(sqrt(a + b*x)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_370():
    f = x**3*(a + b*x)**(sympy.S(1)/3)
    F = -3*a**3*(a + b*x)**(sympy.S(4)/3)/(4*b**4) + 9*a**2*(a + b*x)**(sympy.S(7)/3)/(7*b**4) - 9*a*(a + b*x)**(sympy.S(10)/3)/(10*b**4) + 3*(a + b*x)**(sympy.S(13)/3)/(13*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_371():
    f = x**2*(a + b*x)**(sympy.S(1)/3)
    F = 3*a**2*(a + b*x)**(sympy.S(4)/3)/(4*b**3) - 6*a*(a + b*x)**(sympy.S(7)/3)/(7*b**3) + 3*(a + b*x)**(sympy.S(10)/3)/(10*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_372():
    f = x*(a + b*x)**(sympy.S(1)/3)
    F = -3*a*(a + b*x)**(sympy.S(4)/3)/(4*b**2) + 3*(a + b*x)**(sympy.S(7)/3)/(7*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_373():
    f = (a + b*x)**(sympy.S(1)/3)
    F = 3*(a + b*x)**(sympy.S(4)/3)/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_374():
    f = (a + b*x)**(sympy.S(1)/3)/x
    F = -a**(sympy.S(1)/3)*log(x)/2 + 3*a**(sympy.S(1)/3)*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/2 - sqrt(3)*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3))) + 3*(a + b*x)**(sympy.S(1)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_375():
    f = (a + b*x)**(sympy.S(1)/3)/x**2
    F = -(a + b*x)**(sympy.S(1)/3)/x - b*log(x)/(6*a**(sympy.S(2)/3)) + b*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/(2*a**(sympy.S(2)/3)) - sqrt(3)*b*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_376():
    f = (a + b*x)**(sympy.S(1)/3)/x**3
    F = -(a + b*x)**(sympy.S(1)/3)/(2*x**2) - b*(a + b*x)**(sympy.S(1)/3)/(6*a*x) + b**2*log(x)/(18*a**(sympy.S(5)/3)) - b**2*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/(6*a**(sympy.S(5)/3)) + sqrt(3)*b**2*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_377():
    f = x**3*(a + b*x)**(sympy.S(2)/3)
    F = -3*a**3*(a + b*x)**(sympy.S(5)/3)/(5*b**4) + 9*a**2*(a + b*x)**(sympy.S(8)/3)/(8*b**4) - 9*a*(a + b*x)**(sympy.S(11)/3)/(11*b**4) + 3*(a + b*x)**(sympy.S(14)/3)/(14*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_378():
    f = x**2*(a + b*x)**(sympy.S(2)/3)
    F = 3*a**2*(a + b*x)**(sympy.S(5)/3)/(5*b**3) - 3*a*(a + b*x)**(sympy.S(8)/3)/(4*b**3) + 3*(a + b*x)**(sympy.S(11)/3)/(11*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_379():
    f = x*(a + b*x)**(sympy.S(2)/3)
    F = -3*a*(a + b*x)**(sympy.S(5)/3)/(5*b**2) + 3*(a + b*x)**(sympy.S(8)/3)/(8*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_380():
    f = (a + b*x)**(sympy.S(2)/3)
    F = 3*(a + b*x)**(sympy.S(5)/3)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_381():
    f = (a + b*x)**(sympy.S(2)/3)/x
    F = -a**(sympy.S(2)/3)*log(x)/2 + 3*a**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/2 + sqrt(3)*a**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3))) + 3*(a + b*x)**(sympy.S(2)/3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_382():
    f = (a + b*x)**(sympy.S(2)/3)/x**2
    F = -(a + b*x)**(sympy.S(2)/3)/x - b*log(x)/(3*a**(sympy.S(1)/3)) + b*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/a**(sympy.S(1)/3) + 2*sqrt(3)*b*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_383():
    f = (a + b*x)**(sympy.S(2)/3)/x**3
    F = -(a + b*x)**(sympy.S(2)/3)/(2*x**2) - b*(a + b*x)**(sympy.S(2)/3)/(3*a*x) + b**2*log(x)/(18*a**(sympy.S(4)/3)) - b**2*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/(6*a**(sympy.S(4)/3)) - sqrt(3)*b**2*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_384():
    f = x**3*(a + b*x)**(sympy.S(4)/3)
    F = -3*a**3*(a + b*x)**(sympy.S(7)/3)/(7*b**4) + 9*a**2*(a + b*x)**(sympy.S(10)/3)/(10*b**4) - 9*a*(a + b*x)**(sympy.S(13)/3)/(13*b**4) + 3*(a + b*x)**(sympy.S(16)/3)/(16*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_385():
    f = x**2*(a + b*x)**(sympy.S(4)/3)
    F = 3*a**2*(a + b*x)**(sympy.S(7)/3)/(7*b**3) - 3*a*(a + b*x)**(sympy.S(10)/3)/(5*b**3) + 3*(a + b*x)**(sympy.S(13)/3)/(13*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_386():
    f = x*(a + b*x)**(sympy.S(4)/3)
    F = -3*a*(a + b*x)**(sympy.S(7)/3)/(7*b**2) + 3*(a + b*x)**(sympy.S(10)/3)/(10*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_387():
    f = (a + b*x)**(sympy.S(4)/3)
    F = 3*(a + b*x)**(sympy.S(7)/3)/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_388():
    f = (a + b*x)**(sympy.S(4)/3)/x
    F = -a**(sympy.S(4)/3)*log(x)/2 + 3*a**(sympy.S(4)/3)*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/2 - sqrt(3)*a**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3))) + 3*a*(a + b*x)**(sympy.S(1)/3) + 3*(a + b*x)**(sympy.S(4)/3)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_389():
    f = (a + b*x)**(sympy.S(4)/3)/x**2
    F = -2*a**(sympy.S(1)/3)*b*log(x)/3 + 2*a**(sympy.S(1)/3)*b*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3)) - 4*sqrt(3)*a**(sympy.S(1)/3)*b*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/3 + 4*b*(a + b*x)**(sympy.S(1)/3) - (a + b*x)**(sympy.S(4)/3)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_390():
    f = (a + b*x)**(sympy.S(4)/3)/x**3
    F = -2*b*(a + b*x)**(sympy.S(1)/3)/(3*x) - (a + b*x)**(sympy.S(4)/3)/(2*x**2) - b**2*log(x)/(9*a**(sympy.S(2)/3)) + b**2*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(2)/3)) - 2*sqrt(3)*b**2*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_391():
    f = x**3/(a + b*x)**(sympy.S(1)/3)
    F = -3*a**3*(a + b*x)**(sympy.S(2)/3)/(2*b**4) + 9*a**2*(a + b*x)**(sympy.S(5)/3)/(5*b**4) - 9*a*(a + b*x)**(sympy.S(8)/3)/(8*b**4) + 3*(a + b*x)**(sympy.S(11)/3)/(11*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_392():
    f = x**2/(a + b*x)**(sympy.S(1)/3)
    F = 3*a**2*(a + b*x)**(sympy.S(2)/3)/(2*b**3) - 6*a*(a + b*x)**(sympy.S(5)/3)/(5*b**3) + 3*(a + b*x)**(sympy.S(8)/3)/(8*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_393():
    f = x/(a + b*x)**(sympy.S(1)/3)
    F = -3*a*(a + b*x)**(sympy.S(2)/3)/(2*b**2) + 3*(a + b*x)**(sympy.S(5)/3)/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_394():
    f = (a + b*x)**(sympy.S(-1)/3)
    F = 3*(a + b*x)**(sympy.S(2)/3)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_395():
    f = 1/(x*(a + b*x)**(sympy.S(1)/3))
    F = -log(x)/(2*a**(sympy.S(1)/3)) + 3*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/(2*a**(sympy.S(1)/3)) + sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/a**(sympy.S(1)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_396():
    f = 1/(x**2*(a + b*x)**(sympy.S(1)/3))
    F = -(a + b*x)**(sympy.S(2)/3)/(a*x) + b*log(x)/(6*a**(sympy.S(4)/3)) - b*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/(2*a**(sympy.S(4)/3)) - sqrt(3)*b*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_397():
    f = 1/(x**3*(a + b*x)**(sympy.S(1)/3))
    F = -(a + b*x)**(sympy.S(2)/3)/(2*a*x**2) + 2*b*(a + b*x)**(sympy.S(2)/3)/(3*a**2*x) - b**2*log(x)/(9*a**(sympy.S(7)/3)) + b**2*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(7)/3)) + 2*sqrt(3)*b**2*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_398():
    f = x**3/(-a + b*x)**(sympy.S(1)/3)
    F = 3*a**3*(-a + b*x)**(sympy.S(2)/3)/(2*b**4) + 9*a**2*(-a + b*x)**(sympy.S(5)/3)/(5*b**4) + 9*a*(-a + b*x)**(sympy.S(8)/3)/(8*b**4) + 3*(-a + b*x)**(sympy.S(11)/3)/(11*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_399():
    f = x**2/(-a + b*x)**(sympy.S(1)/3)
    F = 3*a**2*(-a + b*x)**(sympy.S(2)/3)/(2*b**3) + 6*a*(-a + b*x)**(sympy.S(5)/3)/(5*b**3) + 3*(-a + b*x)**(sympy.S(8)/3)/(8*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_400():
    f = x/(-a + b*x)**(sympy.S(1)/3)
    F = 3*a*(-a + b*x)**(sympy.S(2)/3)/(2*b**2) + 3*(-a + b*x)**(sympy.S(5)/3)/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_401():
    f = (-a + b*x)**(sympy.S(-1)/3)
    F = 3*(-a + b*x)**(sympy.S(2)/3)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_402():
    f = 1/(x*(-a + b*x)**(sympy.S(1)/3))
    F = log(x)/(2*a**(sympy.S(1)/3)) - 3*log(a**(sympy.S(1)/3) + (-a + b*x)**(sympy.S(1)/3))/(2*a**(sympy.S(1)/3)) - sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*(-a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/a**(sympy.S(1)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_403():
    f = 1/(x**2*(-a + b*x)**(sympy.S(1)/3))
    F = (-a + b*x)**(sympy.S(2)/3)/(a*x) + b*log(x)/(6*a**(sympy.S(4)/3)) - b*log(a**(sympy.S(1)/3) + (-a + b*x)**(sympy.S(1)/3))/(2*a**(sympy.S(4)/3)) - sqrt(3)*b*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*(-a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_404():
    f = 1/(x**3*(-a + b*x)**(sympy.S(1)/3))
    F = (-a + b*x)**(sympy.S(2)/3)/(2*a*x**2) + 2*b*(-a + b*x)**(sympy.S(2)/3)/(3*a**2*x) + b**2*log(x)/(9*a**(sympy.S(7)/3)) - b**2*log(a**(sympy.S(1)/3) + (-a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(7)/3)) - 2*sqrt(3)*b**2*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*(-a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_405():
    f = x**3/(a + b*x)**(sympy.S(2)/3)
    F = -3*a**3*(a + b*x)**(sympy.S(1)/3)/b**4 + 9*a**2*(a + b*x)**(sympy.S(4)/3)/(4*b**4) - 9*a*(a + b*x)**(sympy.S(7)/3)/(7*b**4) + 3*(a + b*x)**(sympy.S(10)/3)/(10*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_406():
    f = x**2/(a + b*x)**(sympy.S(2)/3)
    F = 3*a**2*(a + b*x)**(sympy.S(1)/3)/b**3 - 3*a*(a + b*x)**(sympy.S(4)/3)/(2*b**3) + 3*(a + b*x)**(sympy.S(7)/3)/(7*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_407():
    f = x/(a + b*x)**(sympy.S(2)/3)
    F = -3*a*(a + b*x)**(sympy.S(1)/3)/b**2 + 3*(a + b*x)**(sympy.S(4)/3)/(4*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_408():
    f = (a + b*x)**(sympy.S(-2)/3)
    F = 3*(a + b*x)**(sympy.S(1)/3)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_409():
    f = 1/(x*(a + b*x)**(sympy.S(2)/3))
    F = -log(x)/(2*a**(sympy.S(2)/3)) + 3*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/(2*a**(sympy.S(2)/3)) - sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/a**(sympy.S(2)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_410():
    f = 1/(x**2*(a + b*x)**(sympy.S(2)/3))
    F = -(a + b*x)**(sympy.S(1)/3)/(a*x) + b*log(x)/(3*a**(sympy.S(5)/3)) - b*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/a**(sympy.S(5)/3) + 2*sqrt(3)*b*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_411():
    f = 1/(x**3*(a + b*x)**(sympy.S(2)/3))
    F = -(a + b*x)**(sympy.S(1)/3)/(2*a*x**2) + 5*b*(a + b*x)**(sympy.S(1)/3)/(6*a**2*x) - 5*b**2*log(x)/(18*a**(sympy.S(8)/3)) + 5*b**2*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/(6*a**(sympy.S(8)/3)) - 5*sqrt(3)*b**2*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_412():
    f = x**3/(a + b*x)**(sympy.S(4)/3)
    F = 3*a**3/(b**4*(a + b*x)**(sympy.S(1)/3)) + 9*a**2*(a + b*x)**(sympy.S(2)/3)/(2*b**4) - 9*a*(a + b*x)**(sympy.S(5)/3)/(5*b**4) + 3*(a + b*x)**(sympy.S(8)/3)/(8*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_413():
    f = x**2/(a + b*x)**(sympy.S(4)/3)
    F = -3*a**2/(b**3*(a + b*x)**(sympy.S(1)/3)) - 3*a*(a + b*x)**(sympy.S(2)/3)/b**3 + 3*(a + b*x)**(sympy.S(5)/3)/(5*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_414():
    f = x/(a + b*x)**(sympy.S(4)/3)
    F = 3*a/(b**2*(a + b*x)**(sympy.S(1)/3)) + 3*(a + b*x)**(sympy.S(2)/3)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_415():
    f = (a + b*x)**(sympy.S(-4)/3)
    F = -3/(b*(a + b*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_416():
    f = 1/(x*(a + b*x)**(sympy.S(4)/3))
    F = 3/(a*(a + b*x)**(sympy.S(1)/3)) - log(x)/(2*a**(sympy.S(4)/3)) + 3*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/(2*a**(sympy.S(4)/3)) + sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/a**(sympy.S(4)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_417():
    f = 1/(x**2*(a + b*x)**(sympy.S(4)/3))
    F = -1/(a*x*(a + b*x)**(sympy.S(1)/3)) - 4*b/(a**2*(a + b*x)**(sympy.S(1)/3)) + 2*b*log(x)/(3*a**(sympy.S(7)/3)) - 2*b*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/a**(sympy.S(7)/3) - 4*sqrt(3)*b*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_418():
    f = 1/(x**3*(a + b*x)**(sympy.S(4)/3))
    F = -1/(2*a*x**2*(a + b*x)**(sympy.S(1)/3)) + 7*b/(6*a**2*x*(a + b*x)**(sympy.S(1)/3)) + 14*b**2/(3*a**3*(a + b*x)**(sympy.S(1)/3)) - 7*b**2*log(x)/(9*a**(sympy.S(10)/3)) + 7*b**2*log(a**(sympy.S(1)/3) - (a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(10)/3)) + 14*sqrt(3)*b**2*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_419():
    f = 1/(x*(a**3 + b**3*x)**(sympy.S(1)/3))
    F = -log(x)/(2*a) + 3*log(a - (a**3 + b**3*x)**(sympy.S(1)/3))/(2*a) + sqrt(3)*atan(sqrt(3)*(a + 2*(a**3 + b**3*x)**(sympy.S(1)/3))/(3*a))/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_420():
    f = 1/(x*(a**3 - b**3*x)**(sympy.S(1)/3))
    F = -log(x)/(2*a) + 3*log(a - (a**3 - b**3*x)**(sympy.S(1)/3))/(2*a) + sqrt(3)*atan(sqrt(3)*(a + 2*(a**3 - b**3*x)**(sympy.S(1)/3))/(3*a))/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_421():
    f = 1/(x*(-a**3 + b**3*x)**(sympy.S(1)/3))
    F = log(x)/(2*a) - 3*log(a + (-a**3 + b**3*x)**(sympy.S(1)/3))/(2*a) - sqrt(3)*atan(sqrt(3)*(a - 2*(-a**3 + b**3*x)**(sympy.S(1)/3))/(3*a))/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_422():
    f = 1/(x*(-a**3 - b**3*x)**(sympy.S(1)/3))
    F = log(x)/(2*a) - 3*log(a + (-a**3 - b**3*x)**(sympy.S(1)/3))/(2*a) - sqrt(3)*atan(sqrt(3)*(a - 2*(-a**3 - b**3*x)**(sympy.S(1)/3))/(3*a))/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_423():
    f = 1/(x*(a**3 + b**3*x)**(sympy.S(2)/3))
    F = -log(x)/(2*a**2) + 3*log(a - (a**3 + b**3*x)**(sympy.S(1)/3))/(2*a**2) - sqrt(3)*atan(sqrt(3)*(a + 2*(a**3 + b**3*x)**(sympy.S(1)/3))/(3*a))/a**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_424():
    f = 1/(x*(a**3 - b**3*x)**(sympy.S(2)/3))
    F = -log(x)/(2*a**2) + 3*log(a - (a**3 - b**3*x)**(sympy.S(1)/3))/(2*a**2) - sqrt(3)*atan(sqrt(3)*(a + 2*(a**3 - b**3*x)**(sympy.S(1)/3))/(3*a))/a**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_425():
    f = 1/(x*(-a**3 + b**3*x)**(sympy.S(2)/3))
    F = -log(x)/(2*a**2) + 3*log(a + (-a**3 + b**3*x)**(sympy.S(1)/3))/(2*a**2) - sqrt(3)*atan(sqrt(3)*(a - 2*(-a**3 + b**3*x)**(sympy.S(1)/3))/(3*a))/a**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_426():
    f = 1/(x*(-a**3 - b**3*x)**(sympy.S(2)/3))
    F = -log(x)/(2*a**2) + 3*log(a + (-a**3 - b**3*x)**(sympy.S(1)/3))/(2*a**2) - sqrt(3)*atan(sqrt(3)*(a - 2*(-a**3 - b**3*x)**(sympy.S(1)/3))/(3*a))/a**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_427():
    f = x**m*(a + b*x)
    F = a*x**(m + 1)/(m + 1) + b*x**(m + 2)/(m + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_428():
    f = x**(sympy.S(5)/2)*(a + b*x)
    F = 2*a*x**(sympy.S(7)/2)/7 + 2*b*x**(sympy.S(9)/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_429():
    f = x**(sympy.S(3)/2)*(a + b*x)
    F = 2*a*x**(sympy.S(5)/2)/5 + 2*b*x**(sympy.S(7)/2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_430():
    f = sqrt(x)*(a + b*x)
    F = 2*a*x**(sympy.S(3)/2)/3 + 2*b*x**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_431():
    f = (a + b*x)/sqrt(x)
    F = 2*a*sqrt(x) + 2*b*x**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_432():
    f = (a + b*x)/x**(sympy.S(3)/2)
    F = -2*a/sqrt(x) + 2*b*sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_433():
    f = (a + b*x)/x**(sympy.S(5)/2)
    F = -2*a/(3*x**(sympy.S(3)/2)) - 2*b/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_434():
    f = x**m*(a + b*x)**2
    F = a**2*x**(m + 1)/(m + 1) + 2*a*b*x**(m + 2)/(m + 2) + b**2*x**(m + 3)/(m + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_435():
    f = x**(sympy.S(5)/2)*(a + b*x)**2
    F = 2*a**2*x**(sympy.S(7)/2)/7 + 4*a*b*x**(sympy.S(9)/2)/9 + 2*b**2*x**(sympy.S(11)/2)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_436():
    f = x**(sympy.S(3)/2)*(a + b*x)**2
    F = 2*a**2*x**(sympy.S(5)/2)/5 + 4*a*b*x**(sympy.S(7)/2)/7 + 2*b**2*x**(sympy.S(9)/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_437():
    f = sqrt(x)*(a + b*x)**2
    F = 2*a**2*x**(sympy.S(3)/2)/3 + 4*a*b*x**(sympy.S(5)/2)/5 + 2*b**2*x**(sympy.S(7)/2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_438():
    f = (a + b*x)**2/sqrt(x)
    F = 2*a**2*sqrt(x) + 4*a*b*x**(sympy.S(3)/2)/3 + 2*b**2*x**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_439():
    f = (a + b*x)**2/x**(sympy.S(3)/2)
    F = -2*a**2/sqrt(x) + 4*a*b*sqrt(x) + 2*b**2*x**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_440():
    f = (a + b*x)**2/x**(sympy.S(5)/2)
    F = -2*a**2/(3*x**(sympy.S(3)/2)) - 4*a*b/sqrt(x) + 2*b**2*sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_441():
    f = x**m*(a + b*x)**3
    F = a**3*x**(m + 1)/(m + 1) + 3*a**2*b*x**(m + 2)/(m + 2) + 3*a*b**2*x**(m + 3)/(m + 3) + b**3*x**(m + 4)/(m + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_442():
    f = x**(sympy.S(5)/2)*(a + b*x)**3
    F = 2*a**3*x**(sympy.S(7)/2)/7 + 2*a**2*b*x**(sympy.S(9)/2)/3 + 6*a*b**2*x**(sympy.S(11)/2)/11 + 2*b**3*x**(sympy.S(13)/2)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_443():
    f = x**(sympy.S(3)/2)*(a + b*x)**3
    F = 2*a**3*x**(sympy.S(5)/2)/5 + 6*a**2*b*x**(sympy.S(7)/2)/7 + 2*a*b**2*x**(sympy.S(9)/2)/3 + 2*b**3*x**(sympy.S(11)/2)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_444():
    f = sqrt(x)*(a + b*x)**3
    F = 2*a**3*x**(sympy.S(3)/2)/3 + 6*a**2*b*x**(sympy.S(5)/2)/5 + 6*a*b**2*x**(sympy.S(7)/2)/7 + 2*b**3*x**(sympy.S(9)/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_445():
    f = (a + b*x)**3/sqrt(x)
    F = 2*a**3*sqrt(x) + 2*a**2*b*x**(sympy.S(3)/2) + 6*a*b**2*x**(sympy.S(5)/2)/5 + 2*b**3*x**(sympy.S(7)/2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_446():
    f = (a + b*x)**3/x**(sympy.S(3)/2)
    F = -2*a**3/sqrt(x) + 6*a**2*b*sqrt(x) + 2*a*b**2*x**(sympy.S(3)/2) + 2*b**3*x**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_447():
    f = (a + b*x)**3/x**(sympy.S(5)/2)
    F = -2*a**3/(3*x**(sympy.S(3)/2)) - 6*a**2*b/sqrt(x) + 6*a*b**2*sqrt(x) + 2*b**3*x**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_448():
    f = x**(sympy.S(5)/2)/(a + b*x)
    F = -2*a**(sympy.S(5)/2)*atan(sqrt(b)*sqrt(x)/sqrt(a))/b**(sympy.S(7)/2) + 2*a**2*sqrt(x)/b**3 - 2*a*x**(sympy.S(3)/2)/(3*b**2) + 2*x**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_449():
    f = x**(sympy.S(3)/2)/(a + b*x)
    F = 2*a**(sympy.S(3)/2)*atan(sqrt(b)*sqrt(x)/sqrt(a))/b**(sympy.S(5)/2) - 2*a*sqrt(x)/b**2 + 2*x**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_450():
    f = sqrt(x)/(a + b*x)
    F = -2*sqrt(a)*atan(sqrt(b)*sqrt(x)/sqrt(a))/b**(sympy.S(3)/2) + 2*sqrt(x)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_451():
    f = 1/(sqrt(x)*(a + b*x))
    F = 2*atan(sqrt(b)*sqrt(x)/sqrt(a))/(sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_452():
    f = 1/(x**(sympy.S(3)/2)*(a + b*x))
    F = -2/(a*sqrt(x)) - 2*sqrt(b)*atan(sqrt(b)*sqrt(x)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_453():
    f = 1/(x**(sympy.S(5)/2)*(a + b*x))
    F = -2/(3*a*x**(sympy.S(3)/2)) + 2*b/(a**2*sqrt(x)) + 2*b**(sympy.S(3)/2)*atan(sqrt(b)*sqrt(x)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_454():
    f = 1/(x**(sympy.S(7)/2)*(a + b*x))
    F = -2/(5*a*x**(sympy.S(5)/2)) + 2*b/(3*a**2*x**(sympy.S(3)/2)) - 2*b**2/(a**3*sqrt(x)) - 2*b**(sympy.S(5)/2)*atan(sqrt(b)*sqrt(x)/sqrt(a))/a**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_455():
    f = x**(sympy.S(5)/2)/(a + b*x)**2
    F = 5*a**(sympy.S(3)/2)*atan(sqrt(b)*sqrt(x)/sqrt(a))/b**(sympy.S(7)/2) - 5*a*sqrt(x)/b**3 - x**(sympy.S(5)/2)/(b*(a + b*x)) + 5*x**(sympy.S(3)/2)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_456():
    f = x**(sympy.S(3)/2)/(a + b*x)**2
    F = -3*sqrt(a)*atan(sqrt(b)*sqrt(x)/sqrt(a))/b**(sympy.S(5)/2) - x**(sympy.S(3)/2)/(b*(a + b*x)) + 3*sqrt(x)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_457():
    f = sqrt(x)/(a + b*x)**2
    F = -sqrt(x)/(b*(a + b*x)) + atan(sqrt(b)*sqrt(x)/sqrt(a))/(sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_458():
    f = 1/(sqrt(x)*(a + b*x)**2)
    F = sqrt(x)/(a*(a + b*x)) + atan(sqrt(b)*sqrt(x)/sqrt(a))/(a**(sympy.S(3)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_459():
    f = 1/(x**(sympy.S(3)/2)*(a + b*x)**2)
    F = 1/(a*sqrt(x)*(a + b*x)) - 3/(a**2*sqrt(x)) - 3*sqrt(b)*atan(sqrt(b)*sqrt(x)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_460():
    f = 1/(x**(sympy.S(5)/2)*(a + b*x)**2)
    F = 1/(a*x**(sympy.S(3)/2)*(a + b*x)) - 5/(3*a**2*x**(sympy.S(3)/2)) + 5*b/(a**3*sqrt(x)) + 5*b**(sympy.S(3)/2)*atan(sqrt(b)*sqrt(x)/sqrt(a))/a**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_461():
    f = x**(sympy.S(7)/2)/(a + b*x)**3
    F = 35*a**(sympy.S(3)/2)*atan(sqrt(b)*sqrt(x)/sqrt(a))/(4*b**(sympy.S(9)/2)) - 35*a*sqrt(x)/(4*b**4) - x**(sympy.S(7)/2)/(2*b*(a + b*x)**2) - 7*x**(sympy.S(5)/2)/(4*b**2*(a + b*x)) + 35*x**(sympy.S(3)/2)/(12*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_462():
    f = x**(sympy.S(5)/2)/(a + b*x)**3
    F = -15*sqrt(a)*atan(sqrt(b)*sqrt(x)/sqrt(a))/(4*b**(sympy.S(7)/2)) - x**(sympy.S(5)/2)/(2*b*(a + b*x)**2) - 5*x**(sympy.S(3)/2)/(4*b**2*(a + b*x)) + 15*sqrt(x)/(4*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_463():
    f = x**(sympy.S(3)/2)/(a + b*x)**3
    F = -x**(sympy.S(3)/2)/(2*b*(a + b*x)**2) - 3*sqrt(x)/(4*b**2*(a + b*x)) + 3*atan(sqrt(b)*sqrt(x)/sqrt(a))/(4*sqrt(a)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_464():
    f = sqrt(x)/(a + b*x)**3
    F = -sqrt(x)/(2*b*(a + b*x)**2) + sqrt(x)/(4*a*b*(a + b*x)) + atan(sqrt(b)*sqrt(x)/sqrt(a))/(4*a**(sympy.S(3)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_465():
    f = 1/(sqrt(x)*(a + b*x)**3)
    F = sqrt(x)/(2*a*(a + b*x)**2) + 3*sqrt(x)/(4*a**2*(a + b*x)) + 3*atan(sqrt(b)*sqrt(x)/sqrt(a))/(4*a**(sympy.S(5)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_466():
    f = 1/(x**(sympy.S(3)/2)*(a + b*x)**3)
    F = 1/(2*a*sqrt(x)*(a + b*x)**2) + 5/(4*a**2*sqrt(x)*(a + b*x)) - 15/(4*a**3*sqrt(x)) - 15*sqrt(b)*atan(sqrt(b)*sqrt(x)/sqrt(a))/(4*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_467():
    f = 1/(x**(sympy.S(5)/2)*(a + b*x)**3)
    F = 1/(2*a*x**(sympy.S(3)/2)*(a + b*x)**2) + 7/(4*a**2*x**(sympy.S(3)/2)*(a + b*x)) - 35/(12*a**3*x**(sympy.S(3)/2)) + 35*b/(4*a**4*sqrt(x)) + 35*b**(sympy.S(3)/2)*atan(sqrt(b)*sqrt(x)/sqrt(a))/(4*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_468():
    f = x**(sympy.S(5)/2)/(-a + b*x)
    F = -2*a**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(a))/b**(sympy.S(7)/2) + 2*a**2*sqrt(x)/b**3 + 2*a*x**(sympy.S(3)/2)/(3*b**2) + 2*x**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_469():
    f = x**(sympy.S(3)/2)/(-a + b*x)
    F = -2*a**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(a))/b**(sympy.S(5)/2) + 2*a*sqrt(x)/b**2 + 2*x**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_470():
    f = sqrt(x)/(-a + b*x)
    F = -2*sqrt(a)*atanh(sqrt(b)*sqrt(x)/sqrt(a))/b**(sympy.S(3)/2) + 2*sqrt(x)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_471():
    f = 1/(sqrt(x)*(-a + b*x))
    F = -2*atanh(sqrt(b)*sqrt(x)/sqrt(a))/(sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_472():
    f = 1/(x**(sympy.S(3)/2)*(-a + b*x))
    F = 2/(a*sqrt(x)) - 2*sqrt(b)*atanh(sqrt(b)*sqrt(x)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_473():
    f = 1/(x**(sympy.S(5)/2)*(-a + b*x))
    F = 2/(3*a*x**(sympy.S(3)/2)) + 2*b/(a**2*sqrt(x)) - 2*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_474():
    f = 1/(x**(sympy.S(7)/2)*(-a + b*x))
    F = 2/(5*a*x**(sympy.S(5)/2)) + 2*b/(3*a**2*x**(sympy.S(3)/2)) + 2*b**2/(a**3*sqrt(x)) - 2*b**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(a))/a**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_475():
    f = x**(sympy.S(5)/2)/(-a + b*x)**2
    F = -5*a**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(a))/b**(sympy.S(7)/2) + 5*a*sqrt(x)/b**3 + x**(sympy.S(5)/2)/(b*(a - b*x)) + 5*x**(sympy.S(3)/2)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_476():
    f = x**(sympy.S(3)/2)/(-a + b*x)**2
    F = -3*sqrt(a)*atanh(sqrt(b)*sqrt(x)/sqrt(a))/b**(sympy.S(5)/2) + x**(sympy.S(3)/2)/(b*(a - b*x)) + 3*sqrt(x)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_477():
    f = sqrt(x)/(-a + b*x)**2
    F = sqrt(x)/(b*(a - b*x)) - atanh(sqrt(b)*sqrt(x)/sqrt(a))/(sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_478():
    f = 1/(sqrt(x)*(-a + b*x)**2)
    F = sqrt(x)/(a*(a - b*x)) + atanh(sqrt(b)*sqrt(x)/sqrt(a))/(a**(sympy.S(3)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_479():
    f = 1/(x**(sympy.S(3)/2)*(-a + b*x)**2)
    F = 1/(a*sqrt(x)*(a - b*x)) - 3/(a**2*sqrt(x)) + 3*sqrt(b)*atanh(sqrt(b)*sqrt(x)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_480():
    f = 1/(x**(sympy.S(5)/2)*(-a + b*x)**2)
    F = 1/(a*x**(sympy.S(3)/2)*(a - b*x)) - 5/(3*a**2*x**(sympy.S(3)/2)) - 5*b/(a**3*sqrt(x)) + 5*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(a))/a**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_481():
    f = x**(sympy.S(7)/2)/(-a + b*x)**3
    F = -35*a**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(a))/(4*b**(sympy.S(9)/2)) + 35*a*sqrt(x)/(4*b**4) - x**(sympy.S(7)/2)/(2*b*(a - b*x)**2) + 7*x**(sympy.S(5)/2)/(4*b**2*(a - b*x)) + 35*x**(sympy.S(3)/2)/(12*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_482():
    f = x**(sympy.S(5)/2)/(-a + b*x)**3
    F = -15*sqrt(a)*atanh(sqrt(b)*sqrt(x)/sqrt(a))/(4*b**(sympy.S(7)/2)) - x**(sympy.S(5)/2)/(2*b*(a - b*x)**2) + 5*x**(sympy.S(3)/2)/(4*b**2*(a - b*x)) + 15*sqrt(x)/(4*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_483():
    f = x**(sympy.S(3)/2)/(-a + b*x)**3
    F = -x**(sympy.S(3)/2)/(2*b*(a - b*x)**2) + 3*sqrt(x)/(4*b**2*(a - b*x)) - 3*atanh(sqrt(b)*sqrt(x)/sqrt(a))/(4*sqrt(a)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_484():
    f = sqrt(x)/(-a + b*x)**3
    F = -sqrt(x)/(2*b*(a - b*x)**2) + sqrt(x)/(4*a*b*(a - b*x)) + atanh(sqrt(b)*sqrt(x)/sqrt(a))/(4*a**(sympy.S(3)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_485():
    f = 1/(sqrt(x)*(-a + b*x)**3)
    F = -sqrt(x)/(2*a*(a - b*x)**2) - 3*sqrt(x)/(4*a**2*(a - b*x)) - 3*atanh(sqrt(b)*sqrt(x)/sqrt(a))/(4*a**(sympy.S(5)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_486():
    f = 1/(x**(sympy.S(3)/2)*(-a + b*x)**3)
    F = -1/(2*a*sqrt(x)*(a - b*x)**2) - 5/(4*a**2*sqrt(x)*(a - b*x)) + 15/(4*a**3*sqrt(x)) - 15*sqrt(b)*atanh(sqrt(b)*sqrt(x)/sqrt(a))/(4*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_487():
    f = 1/(x**(sympy.S(5)/2)*(-a + b*x)**3)
    F = -1/(2*a*x**(sympy.S(3)/2)*(a - b*x)**2) - 7/(4*a**2*x**(sympy.S(3)/2)*(a - b*x)) + 35/(12*a**3*x**(sympy.S(3)/2)) + 35*b/(4*a**4*sqrt(x)) - 35*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(a))/(4*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_488():
    f = x**(sympy.S(5)/2)*sqrt(a + b*x)
    F = -5*a**4*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/(64*b**(sympy.S(7)/2)) + 5*a**3*sqrt(x)*sqrt(a + b*x)/(64*b**3) - 5*a**2*x**(sympy.S(3)/2)*sqrt(a + b*x)/(96*b**2) + a*x**(sympy.S(5)/2)*sqrt(a + b*x)/(24*b) + x**(sympy.S(7)/2)*sqrt(a + b*x)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_489():
    f = x**(sympy.S(3)/2)*sqrt(a + b*x)
    F = a**3*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/(8*b**(sympy.S(5)/2)) - a**2*sqrt(x)*sqrt(a + b*x)/(8*b**2) + a*x**(sympy.S(3)/2)*sqrt(a + b*x)/(12*b) + x**(sympy.S(5)/2)*sqrt(a + b*x)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_490():
    f = sqrt(x)*sqrt(a + b*x)
    F = -a**2*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/(4*b**(sympy.S(3)/2)) + a*sqrt(x)*sqrt(a + b*x)/(4*b) + x**(sympy.S(3)/2)*sqrt(a + b*x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_491():
    f = sqrt(a + b*x)/sqrt(x)
    F = a*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/sqrt(b) + sqrt(x)*sqrt(a + b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_492():
    f = sqrt(a + b*x)/x**(sympy.S(3)/2)
    F = 2*sqrt(b)*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x)) - 2*sqrt(a + b*x)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_493():
    f = sqrt(a + b*x)/x**(sympy.S(5)/2)
    F = -2*(a + b*x)**(sympy.S(3)/2)/(3*a*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_494():
    f = sqrt(a + b*x)/x**(sympy.S(7)/2)
    F = -2*(a + b*x)**(sympy.S(3)/2)/(5*a*x**(sympy.S(5)/2)) + 4*b*(a + b*x)**(sympy.S(3)/2)/(15*a**2*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_495():
    f = sqrt(a + b*x)/x**(sympy.S(9)/2)
    F = -2*(a + b*x)**(sympy.S(3)/2)/(7*a*x**(sympy.S(7)/2)) + 8*b*(a + b*x)**(sympy.S(3)/2)/(35*a**2*x**(sympy.S(5)/2)) - 16*b**2*(a + b*x)**(sympy.S(3)/2)/(105*a**3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_496():
    f = x**(sympy.S(5)/2)*sqrt(a - b*x)
    F = 5*a**4*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/(64*b**(sympy.S(7)/2)) - 5*a**3*sqrt(x)*sqrt(a - b*x)/(64*b**3) - 5*a**2*x**(sympy.S(3)/2)*sqrt(a - b*x)/(96*b**2) - a*x**(sympy.S(5)/2)*sqrt(a - b*x)/(24*b) + x**(sympy.S(7)/2)*sqrt(a - b*x)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_497():
    f = x**(sympy.S(3)/2)*sqrt(a - b*x)
    F = a**3*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/(8*b**(sympy.S(5)/2)) - a**2*sqrt(x)*sqrt(a - b*x)/(8*b**2) - a*x**(sympy.S(3)/2)*sqrt(a - b*x)/(12*b) + x**(sympy.S(5)/2)*sqrt(a - b*x)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_498():
    f = sqrt(x)*sqrt(a - b*x)
    F = a**2*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/(4*b**(sympy.S(3)/2)) - a*sqrt(x)*sqrt(a - b*x)/(4*b) + x**(sympy.S(3)/2)*sqrt(a - b*x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_499():
    f = sqrt(a - b*x)/sqrt(x)
    F = a*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/sqrt(b) + sqrt(x)*sqrt(a - b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_500():
    f = sqrt(a - b*x)/x**(sympy.S(3)/2)
    F = -2*sqrt(b)*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x)) - 2*sqrt(a - b*x)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_501():
    f = sqrt(a - b*x)/x**(sympy.S(5)/2)
    F = -2*(a - b*x)**(sympy.S(3)/2)/(3*a*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_502():
    f = sqrt(a - b*x)/x**(sympy.S(7)/2)
    F = -2*(a - b*x)**(sympy.S(3)/2)/(5*a*x**(sympy.S(5)/2)) - 4*b*(a - b*x)**(sympy.S(3)/2)/(15*a**2*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_503():
    f = sqrt(a - b*x)/x**(sympy.S(9)/2)
    F = -2*(a - b*x)**(sympy.S(3)/2)/(7*a*x**(sympy.S(7)/2)) - 8*b*(a - b*x)**(sympy.S(3)/2)/(35*a**2*x**(sympy.S(5)/2)) - 16*b**2*(a - b*x)**(sympy.S(3)/2)/(105*a**3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_504():
    f = x**(sympy.S(5)/2)*sqrt(b*x + 2)
    F = x**(sympy.S(7)/2)*sqrt(b*x + 2)/4 + x**(sympy.S(5)/2)*sqrt(b*x + 2)/(12*b) - 5*x**(sympy.S(3)/2)*sqrt(b*x + 2)/(24*b**2) + 5*sqrt(x)*sqrt(b*x + 2)/(8*b**3) - 5*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/(4*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_505():
    f = x**(sympy.S(3)/2)*sqrt(b*x + 2)
    F = x**(sympy.S(5)/2)*sqrt(b*x + 2)/3 + x**(sympy.S(3)/2)*sqrt(b*x + 2)/(6*b) - sqrt(x)*sqrt(b*x + 2)/(2*b**2) + asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_506():
    f = sqrt(x)*sqrt(b*x + 2)
    F = x**(sympy.S(3)/2)*sqrt(b*x + 2)/2 + sqrt(x)*sqrt(b*x + 2)/(2*b) - asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_507():
    f = sqrt(b*x + 2)/sqrt(x)
    F = sqrt(x)*sqrt(b*x + 2) + 2*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_508():
    f = sqrt(b*x + 2)/x**(sympy.S(3)/2)
    F = 2*sqrt(b)*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2) - 2*sqrt(b*x + 2)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_509():
    f = sqrt(b*x + 2)/x**(sympy.S(5)/2)
    F = -(b*x + 2)**(sympy.S(3)/2)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_510():
    f = sqrt(b*x + 2)/x**(sympy.S(7)/2)
    F = b*(b*x + 2)**(sympy.S(3)/2)/(15*x**(sympy.S(3)/2)) - (b*x + 2)**(sympy.S(3)/2)/(5*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_511():
    f = sqrt(b*x + 2)/x**(sympy.S(9)/2)
    F = -2*b**2*(b*x + 2)**(sympy.S(3)/2)/(105*x**(sympy.S(3)/2)) + 2*b*(b*x + 2)**(sympy.S(3)/2)/(35*x**(sympy.S(5)/2)) - (b*x + 2)**(sympy.S(3)/2)/(7*x**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_512():
    f = x**(sympy.S(5)/2)*sqrt(-b*x + 2)
    F = x**(sympy.S(7)/2)*sqrt(-b*x + 2)/4 - x**(sympy.S(5)/2)*sqrt(-b*x + 2)/(12*b) - 5*x**(sympy.S(3)/2)*sqrt(-b*x + 2)/(24*b**2) - 5*sqrt(x)*sqrt(-b*x + 2)/(8*b**3) + 5*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/(4*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_513():
    f = x**(sympy.S(3)/2)*sqrt(-b*x + 2)
    F = x**(sympy.S(5)/2)*sqrt(-b*x + 2)/3 - x**(sympy.S(3)/2)*sqrt(-b*x + 2)/(6*b) - sqrt(x)*sqrt(-b*x + 2)/(2*b**2) + asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_514():
    f = sqrt(x)*sqrt(-b*x + 2)
    F = x**(sympy.S(3)/2)*sqrt(-b*x + 2)/2 - sqrt(x)*sqrt(-b*x + 2)/(2*b) + asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_515():
    f = sqrt(-b*x + 2)/sqrt(x)
    F = sqrt(x)*sqrt(-b*x + 2) + 2*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_516():
    f = sqrt(-b*x + 2)/x**(sympy.S(3)/2)
    F = -2*sqrt(b)*asin(sqrt(2)*sqrt(b)*sqrt(x)/2) - 2*sqrt(-b*x + 2)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_517():
    f = sqrt(-b*x + 2)/x**(sympy.S(5)/2)
    F = -(-b*x + 2)**(sympy.S(3)/2)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_518():
    f = sqrt(-b*x + 2)/x**(sympy.S(7)/2)
    F = -b*(-b*x + 2)**(sympy.S(3)/2)/(15*x**(sympy.S(3)/2)) - (-b*x + 2)**(sympy.S(3)/2)/(5*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_519():
    f = sqrt(-b*x + 2)/x**(sympy.S(9)/2)
    F = -2*b**2*(-b*x + 2)**(sympy.S(3)/2)/(105*x**(sympy.S(3)/2)) - 2*b*(-b*x + 2)**(sympy.S(3)/2)/(35*x**(sympy.S(5)/2)) - (-b*x + 2)**(sympy.S(3)/2)/(7*x**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_520():
    f = x**(sympy.S(5)/2)*(a + b*x)**(sympy.S(3)/2)
    F = -3*a**5*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/(128*b**(sympy.S(7)/2)) + 3*a**4*sqrt(x)*sqrt(a + b*x)/(128*b**3) - a**3*x**(sympy.S(3)/2)*sqrt(a + b*x)/(64*b**2) + a**2*x**(sympy.S(5)/2)*sqrt(a + b*x)/(80*b) + 3*a*x**(sympy.S(7)/2)*sqrt(a + b*x)/40 + x**(sympy.S(7)/2)*(a + b*x)**(sympy.S(3)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_521():
    f = x**(sympy.S(3)/2)*(a + b*x)**(sympy.S(3)/2)
    F = 3*a**4*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/(64*b**(sympy.S(5)/2)) - 3*a**3*sqrt(x)*sqrt(a + b*x)/(64*b**2) + a**2*x**(sympy.S(3)/2)*sqrt(a + b*x)/(32*b) + a*x**(sympy.S(5)/2)*sqrt(a + b*x)/8 + x**(sympy.S(5)/2)*(a + b*x)**(sympy.S(3)/2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_522():
    f = sqrt(x)*(a + b*x)**(sympy.S(3)/2)
    F = -a**3*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/(8*b**(sympy.S(3)/2)) + a**2*sqrt(x)*sqrt(a + b*x)/(8*b) + a*x**(sympy.S(3)/2)*sqrt(a + b*x)/4 + x**(sympy.S(3)/2)*(a + b*x)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_523():
    f = (a + b*x)**(sympy.S(3)/2)/sqrt(x)
    F = 3*a**2*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/(4*sqrt(b)) + 3*a*sqrt(x)*sqrt(a + b*x)/4 + sqrt(x)*(a + b*x)**(sympy.S(3)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_524():
    f = (a + b*x)**(sympy.S(3)/2)/x**(sympy.S(3)/2)
    F = 3*a*sqrt(b)*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x)) + 3*b*sqrt(x)*sqrt(a + b*x) - 2*(a + b*x)**(sympy.S(3)/2)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_525():
    f = (a + b*x)**(sympy.S(3)/2)/x**(sympy.S(5)/2)
    F = 2*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x)) - 2*b*sqrt(a + b*x)/sqrt(x) - 2*(a + b*x)**(sympy.S(3)/2)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_526():
    f = x**(sympy.S(5)/2)*(a - b*x)**(sympy.S(3)/2)
    F = 3*a**5*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/(128*b**(sympy.S(7)/2)) - 3*a**4*sqrt(x)*sqrt(a - b*x)/(128*b**3) - a**3*x**(sympy.S(3)/2)*sqrt(a - b*x)/(64*b**2) - a**2*x**(sympy.S(5)/2)*sqrt(a - b*x)/(80*b) + 3*a*x**(sympy.S(7)/2)*sqrt(a - b*x)/40 + x**(sympy.S(7)/2)*(a - b*x)**(sympy.S(3)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_527():
    f = x**(sympy.S(3)/2)*(a - b*x)**(sympy.S(3)/2)
    F = 3*a**4*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/(64*b**(sympy.S(5)/2)) - 3*a**3*sqrt(x)*sqrt(a - b*x)/(64*b**2) - a**2*x**(sympy.S(3)/2)*sqrt(a - b*x)/(32*b) + a*x**(sympy.S(5)/2)*sqrt(a - b*x)/8 + x**(sympy.S(5)/2)*(a - b*x)**(sympy.S(3)/2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_528():
    f = sqrt(x)*(a - b*x)**(sympy.S(3)/2)
    F = a**3*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/(8*b**(sympy.S(3)/2)) - a**2*sqrt(x)*sqrt(a - b*x)/(8*b) + a*x**(sympy.S(3)/2)*sqrt(a - b*x)/4 + x**(sympy.S(3)/2)*(a - b*x)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_529():
    f = (a - b*x)**(sympy.S(3)/2)/sqrt(x)
    F = 3*a**2*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/(4*sqrt(b)) + 3*a*sqrt(x)*sqrt(a - b*x)/4 + sqrt(x)*(a - b*x)**(sympy.S(3)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_530():
    f = (a - b*x)**(sympy.S(3)/2)/x**(sympy.S(3)/2)
    F = -3*a*sqrt(b)*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x)) - 3*b*sqrt(x)*sqrt(a - b*x) - 2*(a - b*x)**(sympy.S(3)/2)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_531():
    f = (a - b*x)**(sympy.S(3)/2)/x**(sympy.S(5)/2)
    F = 2*b**(sympy.S(3)/2)*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x)) + 2*b*sqrt(a - b*x)/sqrt(x) - 2*(a - b*x)**(sympy.S(3)/2)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_532():
    f = x**(sympy.S(5)/2)*(b*x + 2)**(sympy.S(3)/2)
    F = x**(sympy.S(7)/2)*(b*x + 2)**(sympy.S(3)/2)/5 + 3*x**(sympy.S(7)/2)*sqrt(b*x + 2)/20 + x**(sympy.S(5)/2)*sqrt(b*x + 2)/(20*b) - x**(sympy.S(3)/2)*sqrt(b*x + 2)/(8*b**2) + 3*sqrt(x)*sqrt(b*x + 2)/(8*b**3) - 3*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/(4*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_533():
    f = x**(sympy.S(3)/2)*(b*x + 2)**(sympy.S(3)/2)
    F = x**(sympy.S(5)/2)*(b*x + 2)**(sympy.S(3)/2)/4 + x**(sympy.S(5)/2)*sqrt(b*x + 2)/4 + x**(sympy.S(3)/2)*sqrt(b*x + 2)/(8*b) - 3*sqrt(x)*sqrt(b*x + 2)/(8*b**2) + 3*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/(4*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_534():
    f = sqrt(x)*(b*x + 2)**(sympy.S(3)/2)
    F = x**(sympy.S(3)/2)*(b*x + 2)**(sympy.S(3)/2)/3 + x**(sympy.S(3)/2)*sqrt(b*x + 2)/2 + sqrt(x)*sqrt(b*x + 2)/(2*b) - asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_535():
    f = (b*x + 2)**(sympy.S(3)/2)/sqrt(x)
    F = sqrt(x)*(b*x + 2)**(sympy.S(3)/2)/2 + 3*sqrt(x)*sqrt(b*x + 2)/2 + 3*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_536():
    f = (b*x + 2)**(sympy.S(3)/2)/x**(sympy.S(3)/2)
    F = 6*sqrt(b)*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2) + 3*b*sqrt(x)*sqrt(b*x + 2) - 2*(b*x + 2)**(sympy.S(3)/2)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_537():
    f = (b*x + 2)**(sympy.S(3)/2)/x**(sympy.S(5)/2)
    F = 2*b**(sympy.S(3)/2)*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2) - 2*b*sqrt(b*x + 2)/sqrt(x) - 2*(b*x + 2)**(sympy.S(3)/2)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_538():
    f = x**(sympy.S(5)/2)*(-b*x + 2)**(sympy.S(3)/2)
    F = x**(sympy.S(7)/2)*(-b*x + 2)**(sympy.S(3)/2)/5 + 3*x**(sympy.S(7)/2)*sqrt(-b*x + 2)/20 - x**(sympy.S(5)/2)*sqrt(-b*x + 2)/(20*b) - x**(sympy.S(3)/2)*sqrt(-b*x + 2)/(8*b**2) - 3*sqrt(x)*sqrt(-b*x + 2)/(8*b**3) + 3*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/(4*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_539():
    f = x**(sympy.S(3)/2)*(-b*x + 2)**(sympy.S(3)/2)
    F = x**(sympy.S(5)/2)*(-b*x + 2)**(sympy.S(3)/2)/4 + x**(sympy.S(5)/2)*sqrt(-b*x + 2)/4 - x**(sympy.S(3)/2)*sqrt(-b*x + 2)/(8*b) - 3*sqrt(x)*sqrt(-b*x + 2)/(8*b**2) + 3*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/(4*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_540():
    f = sqrt(x)*(-b*x + 2)**(sympy.S(3)/2)
    F = x**(sympy.S(3)/2)*(-b*x + 2)**(sympy.S(3)/2)/3 + x**(sympy.S(3)/2)*sqrt(-b*x + 2)/2 - sqrt(x)*sqrt(-b*x + 2)/(2*b) + asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_541():
    f = (-b*x + 2)**(sympy.S(3)/2)/sqrt(x)
    F = sqrt(x)*(-b*x + 2)**(sympy.S(3)/2)/2 + 3*sqrt(x)*sqrt(-b*x + 2)/2 + 3*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_542():
    f = (-b*x + 2)**(sympy.S(3)/2)/x**(sympy.S(3)/2)
    F = -6*sqrt(b)*asin(sqrt(2)*sqrt(b)*sqrt(x)/2) - 3*b*sqrt(x)*sqrt(-b*x + 2) - 2*(-b*x + 2)**(sympy.S(3)/2)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_543():
    f = (-b*x + 2)**(sympy.S(3)/2)/x**(sympy.S(5)/2)
    F = 2*b**(sympy.S(3)/2)*asin(sqrt(2)*sqrt(b)*sqrt(x)/2) + 2*b*sqrt(-b*x + 2)/sqrt(x) - 2*(-b*x + 2)**(sympy.S(3)/2)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_544():
    f = x**(sympy.S(5)/2)*(a + b*x)**(sympy.S(5)/2)
    F = -5*a**6*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/(512*b**(sympy.S(7)/2)) + 5*a**5*sqrt(x)*sqrt(a + b*x)/(512*b**3) - 5*a**4*x**(sympy.S(3)/2)*sqrt(a + b*x)/(768*b**2) + a**3*x**(sympy.S(5)/2)*sqrt(a + b*x)/(192*b) + a**2*x**(sympy.S(7)/2)*sqrt(a + b*x)/32 + a*x**(sympy.S(7)/2)*(a + b*x)**(sympy.S(3)/2)/12 + x**(sympy.S(7)/2)*(a + b*x)**(sympy.S(5)/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_545():
    f = x**(sympy.S(3)/2)*(a + b*x)**(sympy.S(5)/2)
    F = 3*a**5*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/(128*b**(sympy.S(5)/2)) - 3*a**4*sqrt(x)*sqrt(a + b*x)/(128*b**2) + a**3*x**(sympy.S(3)/2)*sqrt(a + b*x)/(64*b) + a**2*x**(sympy.S(5)/2)*sqrt(a + b*x)/16 + a*x**(sympy.S(5)/2)*(a + b*x)**(sympy.S(3)/2)/8 + x**(sympy.S(5)/2)*(a + b*x)**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_546():
    f = sqrt(x)*(a + b*x)**(sympy.S(5)/2)
    F = -5*a**4*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/(64*b**(sympy.S(3)/2)) + 5*a**3*sqrt(x)*sqrt(a + b*x)/(64*b) + 5*a**2*x**(sympy.S(3)/2)*sqrt(a + b*x)/32 + 5*a*x**(sympy.S(3)/2)*(a + b*x)**(sympy.S(3)/2)/24 + x**(sympy.S(3)/2)*(a + b*x)**(sympy.S(5)/2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_547():
    f = (a + b*x)**(sympy.S(5)/2)/sqrt(x)
    F = 5*a**3*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/(8*sqrt(b)) + 5*a**2*sqrt(x)*sqrt(a + b*x)/8 + 5*a*sqrt(x)*(a + b*x)**(sympy.S(3)/2)/12 + sqrt(x)*(a + b*x)**(sympy.S(5)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_548():
    f = (a + b*x)**(sympy.S(5)/2)/x**(sympy.S(3)/2)
    F = 15*a**2*sqrt(b)*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/4 + 15*a*b*sqrt(x)*sqrt(a + b*x)/4 + 5*b*sqrt(x)*(a + b*x)**(sympy.S(3)/2)/2 - 2*(a + b*x)**(sympy.S(5)/2)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_549():
    f = (a + b*x)**(sympy.S(5)/2)/x**(sympy.S(5)/2)
    F = 5*a*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x)) + 5*b**2*sqrt(x)*sqrt(a + b*x) - 10*b*(a + b*x)**(sympy.S(3)/2)/(3*sqrt(x)) - 2*(a + b*x)**(sympy.S(5)/2)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_550():
    f = x**(sympy.S(5)/2)*(a - b*x)**(sympy.S(5)/2)
    F = 5*a**6*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/(512*b**(sympy.S(7)/2)) - 5*a**5*sqrt(x)*sqrt(a - b*x)/(512*b**3) - 5*a**4*x**(sympy.S(3)/2)*sqrt(a - b*x)/(768*b**2) - a**3*x**(sympy.S(5)/2)*sqrt(a - b*x)/(192*b) + a**2*x**(sympy.S(7)/2)*sqrt(a - b*x)/32 + a*x**(sympy.S(7)/2)*(a - b*x)**(sympy.S(3)/2)/12 + x**(sympy.S(7)/2)*(a - b*x)**(sympy.S(5)/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_551():
    f = x**(sympy.S(3)/2)*(a - b*x)**(sympy.S(5)/2)
    F = 3*a**5*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/(128*b**(sympy.S(5)/2)) - 3*a**4*sqrt(x)*sqrt(a - b*x)/(128*b**2) - a**3*x**(sympy.S(3)/2)*sqrt(a - b*x)/(64*b) + a**2*x**(sympy.S(5)/2)*sqrt(a - b*x)/16 + a*x**(sympy.S(5)/2)*(a - b*x)**(sympy.S(3)/2)/8 + x**(sympy.S(5)/2)*(a - b*x)**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_552():
    f = sqrt(x)*(a - b*x)**(sympy.S(5)/2)
    F = 5*a**4*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/(64*b**(sympy.S(3)/2)) - 5*a**3*sqrt(x)*sqrt(a - b*x)/(64*b) + 5*a**2*x**(sympy.S(3)/2)*sqrt(a - b*x)/32 + 5*a*x**(sympy.S(3)/2)*(a - b*x)**(sympy.S(3)/2)/24 + x**(sympy.S(3)/2)*(a - b*x)**(sympy.S(5)/2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_553():
    f = (a - b*x)**(sympy.S(5)/2)/sqrt(x)
    F = 5*a**3*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/(8*sqrt(b)) + 5*a**2*sqrt(x)*sqrt(a - b*x)/8 + 5*a*sqrt(x)*(a - b*x)**(sympy.S(3)/2)/12 + sqrt(x)*(a - b*x)**(sympy.S(5)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_554():
    f = (a - b*x)**(sympy.S(5)/2)/x**(sympy.S(3)/2)
    F = -15*a**2*sqrt(b)*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/4 - 15*a*b*sqrt(x)*sqrt(a - b*x)/4 - 5*b*sqrt(x)*(a - b*x)**(sympy.S(3)/2)/2 - 2*(a - b*x)**(sympy.S(5)/2)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_555():
    f = (a - b*x)**(sympy.S(5)/2)/x**(sympy.S(5)/2)
    F = 5*a*b**(sympy.S(3)/2)*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x)) + 5*b**2*sqrt(x)*sqrt(a - b*x) + 10*b*(a - b*x)**(sympy.S(3)/2)/(3*sqrt(x)) - 2*(a - b*x)**(sympy.S(5)/2)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_556():
    f = x**(sympy.S(5)/2)*(b*x + 2)**(sympy.S(5)/2)
    F = x**(sympy.S(7)/2)*(b*x + 2)**(sympy.S(5)/2)/6 + x**(sympy.S(7)/2)*(b*x + 2)**(sympy.S(3)/2)/6 + x**(sympy.S(7)/2)*sqrt(b*x + 2)/8 + x**(sympy.S(5)/2)*sqrt(b*x + 2)/(24*b) - 5*x**(sympy.S(3)/2)*sqrt(b*x + 2)/(48*b**2) + 5*sqrt(x)*sqrt(b*x + 2)/(16*b**3) - 5*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/(8*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_557():
    f = x**(sympy.S(3)/2)*(b*x + 2)**(sympy.S(5)/2)
    F = x**(sympy.S(5)/2)*(b*x + 2)**(sympy.S(5)/2)/5 + x**(sympy.S(5)/2)*(b*x + 2)**(sympy.S(3)/2)/4 + x**(sympy.S(5)/2)*sqrt(b*x + 2)/4 + x**(sympy.S(3)/2)*sqrt(b*x + 2)/(8*b) - 3*sqrt(x)*sqrt(b*x + 2)/(8*b**2) + 3*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/(4*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_558():
    f = sqrt(x)*(b*x + 2)**(sympy.S(5)/2)
    F = x**(sympy.S(3)/2)*(b*x + 2)**(sympy.S(5)/2)/4 + 5*x**(sympy.S(3)/2)*(b*x + 2)**(sympy.S(3)/2)/12 + 5*x**(sympy.S(3)/2)*sqrt(b*x + 2)/8 + 5*sqrt(x)*sqrt(b*x + 2)/(8*b) - 5*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/(4*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_559():
    f = (b*x + 2)**(sympy.S(5)/2)/sqrt(x)
    F = sqrt(x)*(b*x + 2)**(sympy.S(5)/2)/3 + 5*sqrt(x)*(b*x + 2)**(sympy.S(3)/2)/6 + 5*sqrt(x)*sqrt(b*x + 2)/2 + 5*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_560():
    f = (b*x + 2)**(sympy.S(5)/2)/x**(sympy.S(3)/2)
    F = 15*sqrt(b)*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2) + 5*b*sqrt(x)*(b*x + 2)**(sympy.S(3)/2)/2 + 15*b*sqrt(x)*sqrt(b*x + 2)/2 - 2*(b*x + 2)**(sympy.S(5)/2)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_561():
    f = (b*x + 2)**(sympy.S(5)/2)/x**(sympy.S(5)/2)
    F = 10*b**(sympy.S(3)/2)*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2) + 5*b**2*sqrt(x)*sqrt(b*x + 2) - 10*b*(b*x + 2)**(sympy.S(3)/2)/(3*sqrt(x)) - 2*(b*x + 2)**(sympy.S(5)/2)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_562():
    f = x**(sympy.S(5)/2)*(-b*x + 2)**(sympy.S(5)/2)
    F = x**(sympy.S(7)/2)*(-b*x + 2)**(sympy.S(5)/2)/6 + x**(sympy.S(7)/2)*(-b*x + 2)**(sympy.S(3)/2)/6 + x**(sympy.S(7)/2)*sqrt(-b*x + 2)/8 - x**(sympy.S(5)/2)*sqrt(-b*x + 2)/(24*b) - 5*x**(sympy.S(3)/2)*sqrt(-b*x + 2)/(48*b**2) - 5*sqrt(x)*sqrt(-b*x + 2)/(16*b**3) + 5*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/(8*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_563():
    f = x**(sympy.S(3)/2)*(-b*x + 2)**(sympy.S(5)/2)
    F = x**(sympy.S(5)/2)*(-b*x + 2)**(sympy.S(5)/2)/5 + x**(sympy.S(5)/2)*(-b*x + 2)**(sympy.S(3)/2)/4 + x**(sympy.S(5)/2)*sqrt(-b*x + 2)/4 - x**(sympy.S(3)/2)*sqrt(-b*x + 2)/(8*b) - 3*sqrt(x)*sqrt(-b*x + 2)/(8*b**2) + 3*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/(4*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_564():
    f = sqrt(x)*(-b*x + 2)**(sympy.S(5)/2)
    F = x**(sympy.S(3)/2)*(-b*x + 2)**(sympy.S(5)/2)/4 + 5*x**(sympy.S(3)/2)*(-b*x + 2)**(sympy.S(3)/2)/12 + 5*x**(sympy.S(3)/2)*sqrt(-b*x + 2)/8 - 5*sqrt(x)*sqrt(-b*x + 2)/(8*b) + 5*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/(4*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_565():
    f = (-b*x + 2)**(sympy.S(5)/2)/sqrt(x)
    F = sqrt(x)*(-b*x + 2)**(sympy.S(5)/2)/3 + 5*sqrt(x)*(-b*x + 2)**(sympy.S(3)/2)/6 + 5*sqrt(x)*sqrt(-b*x + 2)/2 + 5*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_566():
    f = (-b*x + 2)**(sympy.S(5)/2)/x**(sympy.S(3)/2)
    F = -15*sqrt(b)*asin(sqrt(2)*sqrt(b)*sqrt(x)/2) - 5*b*sqrt(x)*(-b*x + 2)**(sympy.S(3)/2)/2 - 15*b*sqrt(x)*sqrt(-b*x + 2)/2 - 2*(-b*x + 2)**(sympy.S(5)/2)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_567():
    f = (-b*x + 2)**(sympy.S(5)/2)/x**(sympy.S(5)/2)
    F = 10*b**(sympy.S(3)/2)*asin(sqrt(2)*sqrt(b)*sqrt(x)/2) + 5*b**2*sqrt(x)*sqrt(-b*x + 2) + 10*b*(-b*x + 2)**(sympy.S(3)/2)/(3*sqrt(x)) - 2*(-b*x + 2)**(sympy.S(5)/2)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_568():
    f = x**(sympy.S(5)/2)/sqrt(a + b*x)
    F = -5*a**3*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/(8*b**(sympy.S(7)/2)) + 5*a**2*sqrt(x)*sqrt(a + b*x)/(8*b**3) - 5*a*x**(sympy.S(3)/2)*sqrt(a + b*x)/(12*b**2) + x**(sympy.S(5)/2)*sqrt(a + b*x)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_569():
    f = x**(sympy.S(3)/2)/sqrt(a + b*x)
    F = 3*a**2*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/(4*b**(sympy.S(5)/2)) - 3*a*sqrt(x)*sqrt(a + b*x)/(4*b**2) + x**(sympy.S(3)/2)*sqrt(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_570():
    f = sqrt(x)/sqrt(a + b*x)
    F = -a*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/b**(sympy.S(3)/2) + sqrt(x)*sqrt(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_571():
    f = 1/(sqrt(x)*sqrt(a + b*x))
    F = 2*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_572():
    f = 1/(x**(sympy.S(3)/2)*sqrt(a + b*x))
    F = -2*sqrt(a + b*x)/(a*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_573():
    f = 1/(x**(sympy.S(5)/2)*sqrt(a + b*x))
    F = -2*sqrt(a + b*x)/(3*a*x**(sympy.S(3)/2)) + 4*b*sqrt(a + b*x)/(3*a**2*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_574():
    f = 1/(x**(sympy.S(7)/2)*sqrt(a + b*x))
    F = -2*sqrt(a + b*x)/(5*a*x**(sympy.S(5)/2)) + 8*b*sqrt(a + b*x)/(15*a**2*x**(sympy.S(3)/2)) - 16*b**2*sqrt(a + b*x)/(15*a**3*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_575():
    f = 1/(x**(sympy.S(9)/2)*sqrt(a + b*x))
    F = -2*sqrt(a + b*x)/(7*a*x**(sympy.S(7)/2)) + 12*b*sqrt(a + b*x)/(35*a**2*x**(sympy.S(5)/2)) - 16*b**2*sqrt(a + b*x)/(35*a**3*x**(sympy.S(3)/2)) + 32*b**3*sqrt(a + b*x)/(35*a**4*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_576():
    f = x**(sympy.S(5)/2)/(a + b*x)**(sympy.S(3)/2)
    F = 15*a**2*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/(4*b**(sympy.S(7)/2)) - 15*a*sqrt(x)*sqrt(a + b*x)/(4*b**3) - 2*x**(sympy.S(5)/2)/(b*sqrt(a + b*x)) + 5*x**(sympy.S(3)/2)*sqrt(a + b*x)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_577():
    f = x**(sympy.S(3)/2)/(a + b*x)**(sympy.S(3)/2)
    F = -3*a*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/b**(sympy.S(5)/2) - 2*x**(sympy.S(3)/2)/(b*sqrt(a + b*x)) + 3*sqrt(x)*sqrt(a + b*x)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_578():
    f = sqrt(x)/(a + b*x)**(sympy.S(3)/2)
    F = -2*sqrt(x)/(b*sqrt(a + b*x)) + 2*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_579():
    f = 1/(sqrt(x)*(a + b*x)**(sympy.S(3)/2))
    F = 2*sqrt(x)/(a*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_580():
    f = 1/(x**(sympy.S(3)/2)*(a + b*x)**(sympy.S(3)/2))
    F = 2/(a*sqrt(x)*sqrt(a + b*x)) - 4*sqrt(a + b*x)/(a**2*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_581():
    f = 1/(x**(sympy.S(5)/2)*(a + b*x)**(sympy.S(3)/2))
    F = 2/(a*x**(sympy.S(3)/2)*sqrt(a + b*x)) - 8*sqrt(a + b*x)/(3*a**2*x**(sympy.S(3)/2)) + 16*b*sqrt(a + b*x)/(3*a**3*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_582():
    f = 1/(x**(sympy.S(7)/2)*(a + b*x)**(sympy.S(3)/2))
    F = 2/(a*x**(sympy.S(5)/2)*sqrt(a + b*x)) - 12*sqrt(a + b*x)/(5*a**2*x**(sympy.S(5)/2)) + 16*b*sqrt(a + b*x)/(5*a**3*x**(sympy.S(3)/2)) - 32*b**2*sqrt(a + b*x)/(5*a**4*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_583():
    f = x**(sympy.S(5)/2)/(a + b*x)**(sympy.S(5)/2)
    F = -5*a*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/b**(sympy.S(7)/2) - 2*x**(sympy.S(5)/2)/(3*b*(a + b*x)**(sympy.S(3)/2)) - 10*x**(sympy.S(3)/2)/(3*b**2*sqrt(a + b*x)) + 5*sqrt(x)*sqrt(a + b*x)/b**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_584():
    f = x**(sympy.S(3)/2)/(a + b*x)**(sympy.S(5)/2)
    F = -2*x**(sympy.S(3)/2)/(3*b*(a + b*x)**(sympy.S(3)/2)) - 2*sqrt(x)/(b**2*sqrt(a + b*x)) + 2*atanh(sqrt(b)*sqrt(x)/sqrt(a + b*x))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_585():
    f = sqrt(x)/(a + b*x)**(sympy.S(5)/2)
    F = 2*x**(sympy.S(3)/2)/(3*a*(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_586():
    f = 1/(sqrt(x)*(a + b*x)**(sympy.S(5)/2))
    F = 2*sqrt(x)/(3*a*(a + b*x)**(sympy.S(3)/2)) + 4*sqrt(x)/(3*a**2*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_587():
    f = 1/(x**(sympy.S(3)/2)*(a + b*x)**(sympy.S(5)/2))
    F = 2/(3*a*sqrt(x)*(a + b*x)**(sympy.S(3)/2)) + 8/(3*a**2*sqrt(x)*sqrt(a + b*x)) - 16*sqrt(a + b*x)/(3*a**3*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_588():
    f = 1/(x**(sympy.S(5)/2)*(a + b*x)**(sympy.S(5)/2))
    F = 2/(3*a*x**(sympy.S(3)/2)*(a + b*x)**(sympy.S(3)/2)) + 4/(a**2*x**(sympy.S(3)/2)*sqrt(a + b*x)) - 16*sqrt(a + b*x)/(3*a**3*x**(sympy.S(3)/2)) + 32*b*sqrt(a + b*x)/(3*a**4*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_589():
    f = x**(sympy.S(5)/2)/sqrt(a - b*x)
    F = 5*a**3*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/(8*b**(sympy.S(7)/2)) - 5*a**2*sqrt(x)*sqrt(a - b*x)/(8*b**3) - 5*a*x**(sympy.S(3)/2)*sqrt(a - b*x)/(12*b**2) - x**(sympy.S(5)/2)*sqrt(a - b*x)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_590():
    f = x**(sympy.S(3)/2)/sqrt(a - b*x)
    F = 3*a**2*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/(4*b**(sympy.S(5)/2)) - 3*a*sqrt(x)*sqrt(a - b*x)/(4*b**2) - x**(sympy.S(3)/2)*sqrt(a - b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_591():
    f = sqrt(x)/sqrt(a - b*x)
    F = a*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/b**(sympy.S(3)/2) - sqrt(x)*sqrt(a - b*x)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_592():
    f = 1/(sqrt(x)*sqrt(a - b*x))
    F = 2*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_593():
    f = 1/(x**(sympy.S(3)/2)*sqrt(a - b*x))
    F = -2*sqrt(a - b*x)/(a*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_594():
    f = 1/(x**(sympy.S(5)/2)*sqrt(a - b*x))
    F = -2*sqrt(a - b*x)/(3*a*x**(sympy.S(3)/2)) - 4*b*sqrt(a - b*x)/(3*a**2*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_595():
    f = x**(sympy.S(5)/2)/(a - b*x)**(sympy.S(3)/2)
    F = -15*a**2*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/(4*b**(sympy.S(7)/2)) + 15*a*sqrt(x)*sqrt(a - b*x)/(4*b**3) + 2*x**(sympy.S(5)/2)/(b*sqrt(a - b*x)) + 5*x**(sympy.S(3)/2)*sqrt(a - b*x)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_596():
    f = x**(sympy.S(3)/2)/(a - b*x)**(sympy.S(3)/2)
    F = -3*a*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/b**(sympy.S(5)/2) + 2*x**(sympy.S(3)/2)/(b*sqrt(a - b*x)) + 3*sqrt(x)*sqrt(a - b*x)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_597():
    f = sqrt(x)/(a - b*x)**(sympy.S(3)/2)
    F = 2*sqrt(x)/(b*sqrt(a - b*x)) - 2*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_598():
    f = 1/(sqrt(x)*(a - b*x)**(sympy.S(3)/2))
    F = 2*sqrt(x)/(a*sqrt(a - b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_599():
    f = 1/(x**(sympy.S(3)/2)*(a - b*x)**(sympy.S(3)/2))
    F = 2/(a*sqrt(x)*sqrt(a - b*x)) - 4*sqrt(a - b*x)/(a**2*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_600():
    f = 1/(x**(sympy.S(5)/2)*(a - b*x)**(sympy.S(3)/2))
    F = 2/(a*x**(sympy.S(3)/2)*sqrt(a - b*x)) - 8*sqrt(a - b*x)/(3*a**2*x**(sympy.S(3)/2)) - 16*b*sqrt(a - b*x)/(3*a**3*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_601():
    f = x**(sympy.S(5)/2)/(a - b*x)**(sympy.S(5)/2)
    F = 5*a*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/b**(sympy.S(7)/2) + 2*x**(sympy.S(5)/2)/(3*b*(a - b*x)**(sympy.S(3)/2)) - 10*x**(sympy.S(3)/2)/(3*b**2*sqrt(a - b*x)) - 5*sqrt(x)*sqrt(a - b*x)/b**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_602():
    f = x**(sympy.S(3)/2)/(a - b*x)**(sympy.S(5)/2)
    F = 2*x**(sympy.S(3)/2)/(3*b*(a - b*x)**(sympy.S(3)/2)) - 2*sqrt(x)/(b**2*sqrt(a - b*x)) + 2*atan(sqrt(b)*sqrt(x)/sqrt(a - b*x))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_603():
    f = sqrt(x)/(a - b*x)**(sympy.S(5)/2)
    F = 2*x**(sympy.S(3)/2)/(3*a*(a - b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_604():
    f = 1/(sqrt(x)*(a - b*x)**(sympy.S(5)/2))
    F = 2*sqrt(x)/(3*a*(a - b*x)**(sympy.S(3)/2)) + 4*sqrt(x)/(3*a**2*sqrt(a - b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_605():
    f = 1/(x**(sympy.S(3)/2)*(a - b*x)**(sympy.S(5)/2))
    F = 2/(3*a*sqrt(x)*(a - b*x)**(sympy.S(3)/2)) + 8/(3*a**2*sqrt(x)*sqrt(a - b*x)) - 16*sqrt(a - b*x)/(3*a**3*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_606():
    f = 1/(x**(sympy.S(5)/2)*(a - b*x)**(sympy.S(5)/2))
    F = 2/(3*a*x**(sympy.S(3)/2)*(a - b*x)**(sympy.S(3)/2)) + 4/(a**2*x**(sympy.S(3)/2)*sqrt(a - b*x)) - 16*sqrt(a - b*x)/(3*a**3*x**(sympy.S(3)/2)) - 32*b*sqrt(a - b*x)/(3*a**4*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_607():
    f = x**(sympy.S(5)/2)/sqrt(b*x + 2)
    F = x**(sympy.S(5)/2)*sqrt(b*x + 2)/(3*b) - 5*x**(sympy.S(3)/2)*sqrt(b*x + 2)/(6*b**2) + 5*sqrt(x)*sqrt(b*x + 2)/(2*b**3) - 5*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_608():
    f = x**(sympy.S(3)/2)/sqrt(b*x + 2)
    F = x**(sympy.S(3)/2)*sqrt(b*x + 2)/(2*b) - 3*sqrt(x)*sqrt(b*x + 2)/(2*b**2) + 3*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_609():
    f = sqrt(x)/sqrt(b*x + 2)
    F = sqrt(x)*sqrt(b*x + 2)/b - 2*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_610():
    f = 1/(sqrt(x)*sqrt(b*x + 2))
    F = 2*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_611():
    f = 1/(x**(sympy.S(3)/2)*sqrt(b*x + 2))
    F = -sqrt(b*x + 2)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_612():
    f = 1/(x**(sympy.S(5)/2)*sqrt(b*x + 2))
    F = b*sqrt(b*x + 2)/(3*sqrt(x)) - sqrt(b*x + 2)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_613():
    f = 1/(x**(sympy.S(7)/2)*sqrt(b*x + 2))
    F = -2*b**2*sqrt(b*x + 2)/(15*sqrt(x)) + 2*b*sqrt(b*x + 2)/(15*x**(sympy.S(3)/2)) - sqrt(b*x + 2)/(5*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_614():
    f = 1/(x**(sympy.S(9)/2)*sqrt(b*x + 2))
    F = 2*b**3*sqrt(b*x + 2)/(35*sqrt(x)) - 2*b**2*sqrt(b*x + 2)/(35*x**(sympy.S(3)/2)) + 3*b*sqrt(b*x + 2)/(35*x**(sympy.S(5)/2)) - sqrt(b*x + 2)/(7*x**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_615():
    f = x**(sympy.S(5)/2)/(b*x + 2)**(sympy.S(3)/2)
    F = -2*x**(sympy.S(5)/2)/(b*sqrt(b*x + 2)) + 5*x**(sympy.S(3)/2)*sqrt(b*x + 2)/(2*b**2) - 15*sqrt(x)*sqrt(b*x + 2)/(2*b**3) + 15*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_616():
    f = x**(sympy.S(3)/2)/(b*x + 2)**(sympy.S(3)/2)
    F = -2*x**(sympy.S(3)/2)/(b*sqrt(b*x + 2)) + 3*sqrt(x)*sqrt(b*x + 2)/b**2 - 6*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_617():
    f = sqrt(x)/(b*x + 2)**(sympy.S(3)/2)
    F = -2*sqrt(x)/(b*sqrt(b*x + 2)) + 2*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_618():
    f = 1/(sqrt(x)*(b*x + 2)**(sympy.S(3)/2))
    F = sqrt(x)/sqrt(b*x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_619():
    f = 1/(x**(sympy.S(3)/2)*(b*x + 2)**(sympy.S(3)/2))
    F = -sqrt(b*x + 2)/sqrt(x) + 1/(sqrt(x)*sqrt(b*x + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_620():
    f = 1/(x**(sympy.S(5)/2)*(b*x + 2)**(sympy.S(3)/2))
    F = 2*b*sqrt(b*x + 2)/(3*sqrt(x)) - 2*sqrt(b*x + 2)/(3*x**(sympy.S(3)/2)) + 1/(x**(sympy.S(3)/2)*sqrt(b*x + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_621():
    f = 1/(x**(sympy.S(7)/2)*(b*x + 2)**(sympy.S(3)/2))
    F = -2*b**2*sqrt(b*x + 2)/(5*sqrt(x)) + 2*b*sqrt(b*x + 2)/(5*x**(sympy.S(3)/2)) - 3*sqrt(b*x + 2)/(5*x**(sympy.S(5)/2)) + 1/(x**(sympy.S(5)/2)*sqrt(b*x + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_622():
    f = x**(sympy.S(5)/2)/(b*x + 2)**(sympy.S(5)/2)
    F = -2*x**(sympy.S(5)/2)/(3*b*(b*x + 2)**(sympy.S(3)/2)) - 10*x**(sympy.S(3)/2)/(3*b**2*sqrt(b*x + 2)) + 5*sqrt(x)*sqrt(b*x + 2)/b**3 - 10*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_623():
    f = x**(sympy.S(3)/2)/(b*x + 2)**(sympy.S(5)/2)
    F = -2*x**(sympy.S(3)/2)/(3*b*(b*x + 2)**(sympy.S(3)/2)) - 2*sqrt(x)/(b**2*sqrt(b*x + 2)) + 2*asinh(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_624():
    f = sqrt(x)/(b*x + 2)**(sympy.S(5)/2)
    F = x**(sympy.S(3)/2)/(3*(b*x + 2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_625():
    f = 1/(sqrt(x)*(b*x + 2)**(sympy.S(5)/2))
    F = sqrt(x)/(3*sqrt(b*x + 2)) + sqrt(x)/(3*(b*x + 2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_626():
    f = 1/(x**(sympy.S(3)/2)*(b*x + 2)**(sympy.S(5)/2))
    F = -2*sqrt(b*x + 2)/(3*sqrt(x)) + 2/(3*sqrt(x)*sqrt(b*x + 2)) + 1/(3*sqrt(x)*(b*x + 2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_627():
    f = 1/(x**(sympy.S(5)/2)*(b*x + 2)**(sympy.S(5)/2))
    F = 2*b*sqrt(b*x + 2)/(3*sqrt(x)) - 2*sqrt(b*x + 2)/(3*x**(sympy.S(3)/2)) + 1/(x**(sympy.S(3)/2)*sqrt(b*x + 2)) + 1/(3*x**(sympy.S(3)/2)*(b*x + 2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_628():
    f = x**(sympy.S(5)/2)/sqrt(-b*x + 2)
    F = -x**(sympy.S(5)/2)*sqrt(-b*x + 2)/(3*b) - 5*x**(sympy.S(3)/2)*sqrt(-b*x + 2)/(6*b**2) - 5*sqrt(x)*sqrt(-b*x + 2)/(2*b**3) + 5*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_629():
    f = x**(sympy.S(3)/2)/sqrt(-b*x + 2)
    F = -x**(sympy.S(3)/2)*sqrt(-b*x + 2)/(2*b) - 3*sqrt(x)*sqrt(-b*x + 2)/(2*b**2) + 3*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_630():
    f = sqrt(x)/sqrt(-b*x + 2)
    F = -sqrt(x)*sqrt(-b*x + 2)/b + 2*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_631():
    f = 1/(sqrt(x)*sqrt(-b*x + 2))
    F = 2*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_632():
    f = 1/(x**(sympy.S(3)/2)*sqrt(-b*x + 2))
    F = -sqrt(-b*x + 2)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_633():
    f = 1/(x**(sympy.S(5)/2)*sqrt(-b*x + 2))
    F = -b*sqrt(-b*x + 2)/(3*sqrt(x)) - sqrt(-b*x + 2)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_634():
    f = x**(sympy.S(5)/2)/(-b*x + 2)**(sympy.S(3)/2)
    F = 2*x**(sympy.S(5)/2)/(b*sqrt(-b*x + 2)) + 5*x**(sympy.S(3)/2)*sqrt(-b*x + 2)/(2*b**2) + 15*sqrt(x)*sqrt(-b*x + 2)/(2*b**3) - 15*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_635():
    f = x**(sympy.S(3)/2)/(-b*x + 2)**(sympy.S(3)/2)
    F = 2*x**(sympy.S(3)/2)/(b*sqrt(-b*x + 2)) + 3*sqrt(x)*sqrt(-b*x + 2)/b**2 - 6*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_636():
    f = sqrt(x)/(-b*x + 2)**(sympy.S(3)/2)
    F = 2*sqrt(x)/(b*sqrt(-b*x + 2)) - 2*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_637():
    f = 1/(sqrt(x)*(-b*x + 2)**(sympy.S(3)/2))
    F = sqrt(x)/sqrt(-b*x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_638():
    f = 1/(x**(sympy.S(3)/2)*(-b*x + 2)**(sympy.S(3)/2))
    F = -sqrt(-b*x + 2)/sqrt(x) + 1/(sqrt(x)*sqrt(-b*x + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_639():
    f = 1/(x**(sympy.S(5)/2)*(-b*x + 2)**(sympy.S(3)/2))
    F = -2*b*sqrt(-b*x + 2)/(3*sqrt(x)) - 2*sqrt(-b*x + 2)/(3*x**(sympy.S(3)/2)) + 1/(x**(sympy.S(3)/2)*sqrt(-b*x + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_640():
    f = x**(sympy.S(5)/2)/(-b*x + 2)**(sympy.S(5)/2)
    F = 2*x**(sympy.S(5)/2)/(3*b*(-b*x + 2)**(sympy.S(3)/2)) - 10*x**(sympy.S(3)/2)/(3*b**2*sqrt(-b*x + 2)) - 5*sqrt(x)*sqrt(-b*x + 2)/b**3 + 10*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_641():
    f = x**(sympy.S(3)/2)/(-b*x + 2)**(sympy.S(5)/2)
    F = 2*x**(sympy.S(3)/2)/(3*b*(-b*x + 2)**(sympy.S(3)/2)) - 2*sqrt(x)/(b**2*sqrt(-b*x + 2)) + 2*asin(sqrt(2)*sqrt(b)*sqrt(x)/2)/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_642():
    f = sqrt(x)/(-b*x + 2)**(sympy.S(5)/2)
    F = x**(sympy.S(3)/2)/(3*(-b*x + 2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_643():
    f = 1/(sqrt(x)*(-b*x + 2)**(sympy.S(5)/2))
    F = sqrt(x)/(3*sqrt(-b*x + 2)) + sqrt(x)/(3*(-b*x + 2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_644():
    f = 1/(x**(sympy.S(3)/2)*(-b*x + 2)**(sympy.S(5)/2))
    F = -2*sqrt(-b*x + 2)/(3*sqrt(x)) + 2/(3*sqrt(x)*sqrt(-b*x + 2)) + 1/(3*sqrt(x)*(-b*x + 2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_645():
    f = 1/(x**(sympy.S(5)/2)*(-b*x + 2)**(sympy.S(5)/2))
    F = -2*b*sqrt(-b*x + 2)/(3*sqrt(x)) - 2*sqrt(-b*x + 2)/(3*x**(sympy.S(3)/2)) + 1/(x**(sympy.S(3)/2)*sqrt(-b*x + 2)) + 1/(3*x**(sympy.S(3)/2)*(-b*x + 2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_646():
    f = sqrt(x)/sqrt(1 - x)
    F = -sqrt(x)*sqrt(1 - x) + asin(2*x - 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_647():
    f = 1/(sqrt(x)*sqrt(1 - x))
    F = asin(2*x - 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_648():
    f = 1/(sqrt(x)*sqrt(-b*x + 1))
    F = 2*asin(sqrt(b)*sqrt(x))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_649():
    f = x**(sympy.S(5)/3)*(a + b*x)
    F = 3*a*x**(sympy.S(8)/3)/8 + 3*b*x**(sympy.S(11)/3)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_650():
    f = x**(sympy.S(4)/3)*(a + b*x)
    F = 3*a*x**(sympy.S(7)/3)/7 + 3*b*x**(sympy.S(10)/3)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_651():
    f = x**(sympy.S(2)/3)*(a + b*x)
    F = 3*a*x**(sympy.S(5)/3)/5 + 3*b*x**(sympy.S(8)/3)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_652():
    f = x**(sympy.S(1)/3)*(a + b*x)
    F = 3*a*x**(sympy.S(4)/3)/4 + 3*b*x**(sympy.S(7)/3)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_653():
    f = (a + b*x)/x**(sympy.S(1)/3)
    F = 3*a*x**(sympy.S(2)/3)/2 + 3*b*x**(sympy.S(5)/3)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_654():
    f = (a + b*x)/x**(sympy.S(2)/3)
    F = 3*a*x**(sympy.S(1)/3) + 3*b*x**(sympy.S(4)/3)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_655():
    f = (a + b*x)/x**(sympy.S(4)/3)
    F = -3*a/x**(sympy.S(1)/3) + 3*b*x**(sympy.S(2)/3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_656():
    f = (a + b*x)/x**(sympy.S(5)/3)
    F = -3*a/(2*x**(sympy.S(2)/3)) + 3*b*x**(sympy.S(1)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_657():
    f = x**(sympy.S(5)/3)*(a + b*x)**2
    F = 3*a**2*x**(sympy.S(8)/3)/8 + 6*a*b*x**(sympy.S(11)/3)/11 + 3*b**2*x**(sympy.S(14)/3)/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_658():
    f = x**(sympy.S(4)/3)*(a + b*x)**2
    F = 3*a**2*x**(sympy.S(7)/3)/7 + 3*a*b*x**(sympy.S(10)/3)/5 + 3*b**2*x**(sympy.S(13)/3)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_659():
    f = x**(sympy.S(2)/3)*(a + b*x)**2
    F = 3*a**2*x**(sympy.S(5)/3)/5 + 3*a*b*x**(sympy.S(8)/3)/4 + 3*b**2*x**(sympy.S(11)/3)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_660():
    f = x**(sympy.S(1)/3)*(a + b*x)**2
    F = 3*a**2*x**(sympy.S(4)/3)/4 + 6*a*b*x**(sympy.S(7)/3)/7 + 3*b**2*x**(sympy.S(10)/3)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_661():
    f = (a + b*x)**2/x**(sympy.S(1)/3)
    F = 3*a**2*x**(sympy.S(2)/3)/2 + 6*a*b*x**(sympy.S(5)/3)/5 + 3*b**2*x**(sympy.S(8)/3)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_662():
    f = (a + b*x)**2/x**(sympy.S(2)/3)
    F = 3*a**2*x**(sympy.S(1)/3) + 3*a*b*x**(sympy.S(4)/3)/2 + 3*b**2*x**(sympy.S(7)/3)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_663():
    f = (a + b*x)**2/x**(sympy.S(4)/3)
    F = -3*a**2/x**(sympy.S(1)/3) + 3*a*b*x**(sympy.S(2)/3) + 3*b**2*x**(sympy.S(5)/3)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_664():
    f = (a + b*x)**2/x**(sympy.S(5)/3)
    F = -3*a**2/(2*x**(sympy.S(2)/3)) + 6*a*b*x**(sympy.S(1)/3) + 3*b**2*x**(sympy.S(4)/3)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_665():
    f = x**(sympy.S(5)/3)*(a + b*x)**3
    F = 3*a**3*x**(sympy.S(8)/3)/8 + 9*a**2*b*x**(sympy.S(11)/3)/11 + 9*a*b**2*x**(sympy.S(14)/3)/14 + 3*b**3*x**(sympy.S(17)/3)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_666():
    f = x**(sympy.S(4)/3)*(a + b*x)**3
    F = 3*a**3*x**(sympy.S(7)/3)/7 + 9*a**2*b*x**(sympy.S(10)/3)/10 + 9*a*b**2*x**(sympy.S(13)/3)/13 + 3*b**3*x**(sympy.S(16)/3)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_667():
    f = x**(sympy.S(2)/3)*(a + b*x)**3
    F = 3*a**3*x**(sympy.S(5)/3)/5 + 9*a**2*b*x**(sympy.S(8)/3)/8 + 9*a*b**2*x**(sympy.S(11)/3)/11 + 3*b**3*x**(sympy.S(14)/3)/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_668():
    f = x**(sympy.S(1)/3)*(a + b*x)**3
    F = 3*a**3*x**(sympy.S(4)/3)/4 + 9*a**2*b*x**(sympy.S(7)/3)/7 + 9*a*b**2*x**(sympy.S(10)/3)/10 + 3*b**3*x**(sympy.S(13)/3)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_669():
    f = (a + b*x)**3/x**(sympy.S(1)/3)
    F = 3*a**3*x**(sympy.S(2)/3)/2 + 9*a**2*b*x**(sympy.S(5)/3)/5 + 9*a*b**2*x**(sympy.S(8)/3)/8 + 3*b**3*x**(sympy.S(11)/3)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_670():
    f = (a + b*x)**3/x**(sympy.S(2)/3)
    F = 3*a**3*x**(sympy.S(1)/3) + 9*a**2*b*x**(sympy.S(4)/3)/4 + 9*a*b**2*x**(sympy.S(7)/3)/7 + 3*b**3*x**(sympy.S(10)/3)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_671():
    f = (a + b*x)**3/x**(sympy.S(4)/3)
    F = -3*a**3/x**(sympy.S(1)/3) + 9*a**2*b*x**(sympy.S(2)/3)/2 + 9*a*b**2*x**(sympy.S(5)/3)/5 + 3*b**3*x**(sympy.S(8)/3)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_672():
    f = (a + b*x)**3/x**(sympy.S(5)/3)
    F = -3*a**3/(2*x**(sympy.S(2)/3)) + 9*a**2*b*x**(sympy.S(1)/3) + 9*a*b**2*x**(sympy.S(4)/3)/4 + 3*b**3*x**(sympy.S(7)/3)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_673():
    f = x**(sympy.S(5)/3)/(a + b*x)
    F = -3*a**(sympy.S(5)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(2*b**(sympy.S(8)/3)) + a**(sympy.S(5)/3)*log(a + b*x)/(2*b**(sympy.S(8)/3)) - sqrt(3)*a**(sympy.S(5)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/b**(sympy.S(8)/3) - 3*a*x**(sympy.S(2)/3)/(2*b**2) + 3*x**(sympy.S(5)/3)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_674():
    f = x**(sympy.S(4)/3)/(a + b*x)
    F = 3*a**(sympy.S(4)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(2*b**(sympy.S(7)/3)) - a**(sympy.S(4)/3)*log(a + b*x)/(2*b**(sympy.S(7)/3)) - sqrt(3)*a**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/b**(sympy.S(7)/3) - 3*a*x**(sympy.S(1)/3)/b**2 + 3*x**(sympy.S(4)/3)/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_675():
    f = x**(sympy.S(2)/3)/(a + b*x)
    F = 3*a**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(2*b**(sympy.S(5)/3)) - a**(sympy.S(2)/3)*log(a + b*x)/(2*b**(sympy.S(5)/3)) + sqrt(3)*a**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/b**(sympy.S(5)/3) + 3*x**(sympy.S(2)/3)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_676():
    f = x**(sympy.S(1)/3)/(a + b*x)
    F = -3*a**(sympy.S(1)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(2*b**(sympy.S(4)/3)) + a**(sympy.S(1)/3)*log(a + b*x)/(2*b**(sympy.S(4)/3)) + sqrt(3)*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/b**(sympy.S(4)/3) + 3*x**(sympy.S(1)/3)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_677():
    f = 1/(x**(sympy.S(1)/3)*(a + b*x))
    F = -3*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(2*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)) + log(a + b*x)/(2*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)) - sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(a**(sympy.S(1)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_678():
    f = 1/(x**(sympy.S(2)/3)*(a + b*x))
    F = 3*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(2*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)) - log(a + b*x)/(2*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)) - sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(a**(sympy.S(2)/3)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_679():
    f = 1/(x**(sympy.S(4)/3)*(a + b*x))
    F = -3/(a*x**(sympy.S(1)/3)) + 3*b**(sympy.S(1)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(2*a**(sympy.S(4)/3)) - b**(sympy.S(1)/3)*log(a + b*x)/(2*a**(sympy.S(4)/3)) + sqrt(3)*b**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/a**(sympy.S(4)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_680():
    f = 1/(x**(sympy.S(5)/3)*(a + b*x))
    F = -3/(2*a*x**(sympy.S(2)/3)) - 3*b**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(2*a**(sympy.S(5)/3)) + b**(sympy.S(2)/3)*log(a + b*x)/(2*a**(sympy.S(5)/3)) + sqrt(3)*b**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/a**(sympy.S(5)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_681():
    f = x**(sympy.S(5)/3)/(a + b*x)**2
    F = 5*a**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(2*b**(sympy.S(8)/3)) - 5*a**(sympy.S(2)/3)*log(a + b*x)/(6*b**(sympy.S(8)/3)) + 5*sqrt(3)*a**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(8)/3)) - x**(sympy.S(5)/3)/(b*(a + b*x)) + 5*x**(sympy.S(2)/3)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_682():
    f = x**(sympy.S(4)/3)/(a + b*x)**2
    F = -2*a**(sympy.S(1)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/b**(sympy.S(7)/3) + 2*a**(sympy.S(1)/3)*log(a + b*x)/(3*b**(sympy.S(7)/3)) + 4*sqrt(3)*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(7)/3)) - x**(sympy.S(4)/3)/(b*(a + b*x)) + 4*x**(sympy.S(1)/3)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_683():
    f = x**(sympy.S(2)/3)/(a + b*x)**2
    F = -x**(sympy.S(2)/3)/(b*(a + b*x)) - log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*b**(sympy.S(5)/3)) + log(a + b*x)/(3*a**(sympy.S(1)/3)*b**(sympy.S(5)/3)) - 2*sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_684():
    f = x**(sympy.S(1)/3)/(a + b*x)**2
    F = -x**(sympy.S(1)/3)/(b*(a + b*x)) + log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(2*a**(sympy.S(2)/3)*b**(sympy.S(4)/3)) - log(a + b*x)/(6*a**(sympy.S(2)/3)*b**(sympy.S(4)/3)) - sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_685():
    f = 1/(x**(sympy.S(1)/3)*(a + b*x)**2)
    F = x**(sympy.S(2)/3)/(a*(a + b*x)) - log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(2*a**(sympy.S(4)/3)*b**(sympy.S(2)/3)) + log(a + b*x)/(6*a**(sympy.S(4)/3)*b**(sympy.S(2)/3)) - sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(4)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_686():
    f = 1/(x**(sympy.S(2)/3)*(a + b*x)**2)
    F = x**(sympy.S(1)/3)/(a*(a + b*x)) + log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(a**(sympy.S(5)/3)*b**(sympy.S(1)/3)) - log(a + b*x)/(3*a**(sympy.S(5)/3)*b**(sympy.S(1)/3)) - 2*sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/3)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_687():
    f = 1/(x**(sympy.S(4)/3)*(a + b*x)**2)
    F = 1/(a*x**(sympy.S(1)/3)*(a + b*x)) - 4/(a**2*x**(sympy.S(1)/3)) + 2*b**(sympy.S(1)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/a**(sympy.S(7)/3) - 2*b**(sympy.S(1)/3)*log(a + b*x)/(3*a**(sympy.S(7)/3)) + 4*sqrt(3)*b**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_688():
    f = 1/(x**(sympy.S(5)/3)*(a + b*x)**2)
    F = 1/(a*x**(sympy.S(2)/3)*(a + b*x)) - 5/(2*a**2*x**(sympy.S(2)/3)) - 5*b**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(2*a**(sympy.S(8)/3)) + 5*b**(sympy.S(2)/3)*log(a + b*x)/(6*a**(sympy.S(8)/3)) + 5*sqrt(3)*b**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_689():
    f = x**(sympy.S(5)/3)/(a + b*x)**3
    F = -x**(sympy.S(5)/3)/(2*b*(a + b*x)**2) - 5*x**(sympy.S(2)/3)/(6*b**2*(a + b*x)) - 5*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(6*a**(sympy.S(1)/3)*b**(sympy.S(8)/3)) + 5*log(a + b*x)/(18*a**(sympy.S(1)/3)*b**(sympy.S(8)/3)) - 5*sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(1)/3)*b**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_690():
    f = x**(sympy.S(4)/3)/(a + b*x)**3
    F = -x**(sympy.S(4)/3)/(2*b*(a + b*x)**2) - 2*x**(sympy.S(1)/3)/(3*b**2*(a + b*x)) + log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(2)/3)*b**(sympy.S(7)/3)) - log(a + b*x)/(9*a**(sympy.S(2)/3)*b**(sympy.S(7)/3)) - 2*sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(2)/3)*b**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_691():
    f = x**(sympy.S(2)/3)/(a + b*x)**3
    F = -x**(sympy.S(2)/3)/(2*b*(a + b*x)**2) + x**(sympy.S(2)/3)/(3*a*b*(a + b*x)) - log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(6*a**(sympy.S(4)/3)*b**(sympy.S(5)/3)) + log(a + b*x)/(18*a**(sympy.S(4)/3)*b**(sympy.S(5)/3)) - sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(4)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_692():
    f = x**(sympy.S(1)/3)/(a + b*x)**3
    F = -x**(sympy.S(1)/3)/(2*b*(a + b*x)**2) + x**(sympy.S(1)/3)/(6*a*b*(a + b*x)) + log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(6*a**(sympy.S(5)/3)*b**(sympy.S(4)/3)) - log(a + b*x)/(18*a**(sympy.S(5)/3)*b**(sympy.S(4)/3)) - sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_693():
    f = 1/(x**(sympy.S(1)/3)*(a + b*x)**3)
    F = x**(sympy.S(2)/3)/(2*a*(a + b*x)**2) + 2*x**(sympy.S(2)/3)/(3*a**2*(a + b*x)) - log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(7)/3)*b**(sympy.S(2)/3)) + log(a + b*x)/(9*a**(sympy.S(7)/3)*b**(sympy.S(2)/3)) - 2*sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(7)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_694():
    f = 1/(x**(sympy.S(2)/3)*(a + b*x)**3)
    F = x**(sympy.S(1)/3)/(2*a*(a + b*x)**2) + 5*x**(sympy.S(1)/3)/(6*a**2*(a + b*x)) + 5*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(6*a**(sympy.S(8)/3)*b**(sympy.S(1)/3)) - 5*log(a + b*x)/(18*a**(sympy.S(8)/3)*b**(sympy.S(1)/3)) - 5*sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(8)/3)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_695():
    f = 1/(x**(sympy.S(4)/3)*(a + b*x)**3)
    F = 1/(2*a*x**(sympy.S(1)/3)*(a + b*x)**2) + 7/(6*a**2*x**(sympy.S(1)/3)*(a + b*x)) - 14/(3*a**3*x**(sympy.S(1)/3)) + 7*b**(sympy.S(1)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(10)/3)) - 7*b**(sympy.S(1)/3)*log(a + b*x)/(9*a**(sympy.S(10)/3)) + 14*sqrt(3)*b**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_696():
    f = 1/(x**(sympy.S(5)/3)*(a + b*x)**3)
    F = 1/(2*a*x**(sympy.S(2)/3)*(a + b*x)**2) + 4/(3*a**2*x**(sympy.S(2)/3)*(a + b*x)) - 10/(3*a**3*x**(sympy.S(2)/3)) - 10*b**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(11)/3)) + 10*b**(sympy.S(2)/3)*log(a + b*x)/(9*a**(sympy.S(11)/3)) + 20*sqrt(3)*b**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(11)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_697():
    f = (1 - x)**(sympy.S(1)/4)/(x + 1)
    F = 4*(1 - x)**(sympy.S(1)/4) - 2*2**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*(1 - x)**(sympy.S(1)/4)/2) - 2*2**(sympy.S(1)/4)*atanh(2**(sympy.S(3)/4)*(1 - x)**(sympy.S(1)/4)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_698():
    f = x**m*(a + b*x)**10
    F = a**10*x**(m + 1)/(m + 1) + 10*a**9*b*x**(m + 2)/(m + 2) + 45*a**8*b**2*x**(m + 3)/(m + 3) + 120*a**7*b**3*x**(m + 4)/(m + 4) + 210*a**6*b**4*x**(m + 5)/(m + 5) + 252*a**5*b**5*x**(m + 6)/(m + 6) + 210*a**4*b**6*x**(m + 7)/(m + 7) + 120*a**3*b**7*x**(m + 8)/(m + 8) + 45*a**2*b**8*x**(m + 9)/(m + 9) + 10*a*b**9*x**(m + 10)/(m + 10) + b**10*x**(m + 11)/(m + 11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_699():
    f = x**m*(a + b*x)**7
    F = a**7*x**(m + 1)/(m + 1) + 7*a**6*b*x**(m + 2)/(m + 2) + 21*a**5*b**2*x**(m + 3)/(m + 3) + 35*a**4*b**3*x**(m + 4)/(m + 4) + 35*a**3*b**4*x**(m + 5)/(m + 5) + 21*a**2*b**5*x**(m + 6)/(m + 6) + 7*a*b**6*x**(m + 7)/(m + 7) + b**7*x**(m + 8)/(m + 8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_700():
    f = x**m*(a + b*x)**3
    F = a**3*x**(m + 1)/(m + 1) + 3*a**2*b*x**(m + 2)/(m + 2) + 3*a*b**2*x**(m + 3)/(m + 3) + b**3*x**(m + 4)/(m + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_701():
    f = x**m*(a + b*x)**2
    F = a**2*x**(m + 1)/(m + 1) + 2*a*b*x**(m + 2)/(m + 2) + b**2*x**(m + 3)/(m + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_702():
    f = x**m*(a + b*x)
    F = a*x**(m + 1)/(m + 1) + b*x**(m + 2)/(m + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_703():
    f = x**m/(a + b*x)
    F = x**(m + 1)*hyper((1, m + 1), (m + 2,), -b*x/a)/(a*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_704():
    f = x**m/(a + b*x)**2
    F = x**(m + 1)*hyper((2, m + 1), (m + 2,), -b*x/a)/(a**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_705():
    f = x**m/(a + b*x)**3
    F = x**(m + 1)*hyper((3, m + 1), (m + 2,), -b*x/a)/(a**3*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_706():
    f = x**m*(a + b*x)**(sympy.S(5)/2)
    F = 2*x**m*(a + b*x)**(sympy.S(7)/2)*hyper((sympy.S(7)/2, -m), (sympy.S(9)/2,), 1 + b*x/a)/(7*b*(-b*x/a)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_707():
    f = x**m*(a + b*x)**(sympy.S(3)/2)
    F = 2*x**m*(a + b*x)**(sympy.S(5)/2)*hyper((sympy.S(5)/2, -m), (sympy.S(7)/2,), 1 + b*x/a)/(5*b*(-b*x/a)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_708():
    f = x**m*sqrt(a + b*x)
    F = 2*x**m*(a + b*x)**(sympy.S(3)/2)*hyper((sympy.S(3)/2, -m), (sympy.S(5)/2,), 1 + b*x/a)/(3*b*(-b*x/a)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_709():
    f = x**m/sqrt(a + b*x)
    F = 2*x**m*sqrt(a + b*x)*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), 1 + b*x/a)/(b*(-b*x/a)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_710():
    f = x**m/(a + b*x)**(sympy.S(3)/2)
    F = -2*x**m*hyper((sympy.S(-1)/2, -m), (sympy.S.Half,), 1 + b*x/a)/(b*(-b*x/a)**m*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_711():
    f = x**m/(a + b*x)**(sympy.S(5)/2)
    F = -2*x**m*hyper((sympy.S(-3)/2, -m), (sympy.S(-1)/2,), 1 + b*x/a)/(3*b*(-b*x/a)**m*(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_712():
    f = x**(m + 2)/sqrt(a + b*x)
    F = 2*a**2*x**m*sqrt(a + b*x)*hyper((sympy.S.Half, -m - 2), (sympy.S(3)/2,), 1 + b*x/a)/(b**3*(-b*x/a)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_713():
    f = x**(m + 1)/sqrt(a + b*x)
    F = -2*a*x**m*sqrt(a + b*x)*hyper((sympy.S.Half, -m - 1), (sympy.S(3)/2,), 1 + b*x/a)/(b**2*(-b*x/a)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_714():
    f = x**m/sqrt(a + b*x)
    F = 2*x**m*sqrt(a + b*x)*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), 1 + b*x/a)/(b*(-b*x/a)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_715():
    f = x**(m - 1)/sqrt(a + b*x)
    F = -2*x**m*sqrt(a + b*x)*hyper((sympy.S.Half, 1 - m), (sympy.S(3)/2,), 1 + b*x/a)/(a*(-b*x/a)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_716():
    f = x**(m - 2)/sqrt(a + b*x)
    F = 2*b*x**m*sqrt(a + b*x)*hyper((sympy.S.Half, 2 - m), (sympy.S(3)/2,), 1 + b*x/a)/(a**2*(-b*x/a)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_717():
    f = x**(m - 3)/sqrt(a + b*x)
    F = -2*b**2*x**m*sqrt(a + b*x)*hyper((sympy.S.Half, 3 - m), (sympy.S(3)/2,), 1 + b*x/a)/(a**3*(-b*x/a)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_718():
    f = x**m/sqrt(3*x + 2)
    F = sqrt(2)*x**(m + 1)*hyper((sympy.S.Half, m + 1), (m + 2,), -3*x/2)/(2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_719():
    f = x**m/sqrt(2 - 3*x)
    F = sqrt(2)*x**(m + 1)*hyper((sympy.S.Half, m + 1), (m + 2,), 3*x/2)/(2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_720():
    f = x**m/sqrt(3*x - 2)
    F = (sympy.S(3)/2)**(-m - 1)*sqrt(3*x - 2)*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), 1 - 3*x/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_721():
    f = x**m/sqrt(-3*x - 2)
    F = -2**(m + 1)*3**(-m - 1)*x**m*sqrt(-3*x - 2)*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), 3*x/2 + 1)/(-x)**m
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_722():
    f = (-x)**m/sqrt(a + b*x)
    F = 2*(-x)**m*sqrt(a + b*x)*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), 1 + b*x/a)/(b*(-b*x/a)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_723():
    f = (-x)**m/sqrt(3*x + 2)
    F = -sqrt(2)*(-x)**(m + 1)*hyper((sympy.S.Half, m + 1), (m + 2,), -3*x/2)/(2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_724():
    f = (-x)**m/sqrt(2 - 3*x)
    F = -sqrt(2)*(-x)**(m + 1)*hyper((sympy.S.Half, m + 1), (m + 2,), 3*x/2)/(2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_725():
    f = (-x)**m/sqrt(3*x - 2)
    F = 2**(m + 1)*3**(-m - 1)*(-x)**m*sqrt(3*x - 2)*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), 1 - 3*x/2)/x**m
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_726():
    f = (-x)**m/sqrt(-3*x - 2)
    F = -(sympy.S(3)/2)**(-m - 1)*sqrt(-3*x - 2)*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), 3*x/2 + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_727():
    f = x**n/sqrt(1 - x)
    F = -2*sqrt(1 - x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), 1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_728():
    f = x**n/sqrt(-a*x + a)
    F = -2*sqrt(-a*x + a)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), 1 - x)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_729():
    f = x**m*(a + b*x)**n
    F = x**(m + 1)*(a + b*x)**n*hyper((-n, m + 1), (m + 2,), -b*x/a)/((1 + b*x/a)**n*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_730():
    f = (c*x)**m*(a + b*x)**n
    F = (c*x)**(m + 1)*(a + b*x)**n*hyper((-n, m + 1), (m + 2,), -b*x/a)/(c*(1 + b*x/a)**n*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_731():
    f = x**3*(a + b*x)**n
    F = -a**3*(a + b*x)**(n + 1)/(b**4*(n + 1)) + 3*a**2*(a + b*x)**(n + 2)/(b**4*(n + 2)) - 3*a*(a + b*x)**(n + 3)/(b**4*(n + 3)) + (a + b*x)**(n + 4)/(b**4*(n + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_732():
    f = x**2*(a + b*x)**n
    F = a**2*(a + b*x)**(n + 1)/(b**3*(n + 1)) - 2*a*(a + b*x)**(n + 2)/(b**3*(n + 2)) + (a + b*x)**(n + 3)/(b**3*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_733():
    f = x*(a + b*x)**n
    F = -a*(a + b*x)**(n + 1)/(b**2*(n + 1)) + (a + b*x)**(n + 2)/(b**2*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_734():
    f = (a + b*x)**n
    F = (a + b*x)**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_735():
    f = (a + b*x)**n/x
    F = -(a + b*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*x/a)/(a*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_736():
    f = (a + b*x)**n/x**2
    F = b*(a + b*x)**(n + 1)*hyper((2, n + 1), (n + 2,), 1 + b*x/a)/(a**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_737():
    f = (a + b*x)**n/x**3
    F = -b**2*(a + b*x)**(n + 1)*hyper((3, n + 1), (n + 2,), 1 + b*x/a)/(a**3*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_738():
    f = x**(n - 4)/(a + b*x)**n
    F = -x**(n - 3)*(a + b*x)**(1 - n)/(a*(3 - n)) + 2*b*x**(n - 2)*(a + b*x)**(1 - n)/(a**2*(2 - n)*(3 - n)) - 2*b**2*x**(n - 1)*(a + b*x)**(1 - n)/(a**3*(1 - n)*(2 - n)*(3 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_739():
    f = x**(n - 3)/(a + b*x)**n
    F = -x**(n - 2)*(a + b*x)**(1 - n)/(a*(2 - n)) + b*x**(n - 1)*(a + b*x)**(1 - n)/(a**2*(1 - n)*(2 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_740():
    f = x**(n - 2)/(a + b*x)**n
    F = -x**(n - 1)*(a + b*x)**(1 - n)/(a*(1 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_741():
    f = x**(n - 1)/(a + b*x)**n
    F = x**n*(1 + b*x/a)**n*hyper((n, n), (n + 1,), -b*x/a)/(n*(a + b*x)**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_742():
    f = x**n/(a + b*x)**n
    F = x**(n + 1)*(1 + b*x/a)**n*hyper((n, n + 1), (n + 2,), -b*x/a)/((a + b*x)**n*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_743():
    f = x**(n + 1)/(a + b*x)**n
    F = x**(n + 2)*(1 + b*x/a)**n*hyper((n, n + 2), (n + 3,), -b*x/a)/((a + b*x)**n*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_744():
    f = x**(sympy.S(3)/2)*(a + b*x)**n
    F = 2*x**(sympy.S(5)/2)*(a + b*x)**n*hyper((sympy.S(5)/2, -n), (sympy.S(7)/2,), -b*x/a)/(5*(1 + b*x/a)**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_745():
    f = sqrt(x)*(a + b*x)**n
    F = 2*x**(sympy.S(3)/2)*(a + b*x)**n*hyper((sympy.S(3)/2, -n), (sympy.S(5)/2,), -b*x/a)/(3*(1 + b*x/a)**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_746():
    f = (a + b*x)**n/sqrt(x)
    F = 2*sqrt(x)*(a + b*x)**n*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), -b*x/a)/(1 + b*x/a)**n
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_747():
    f = (a + b*x)**n/x**(sympy.S(3)/2)
    F = -2*(a + b*x)**n*hyper((sympy.S(-1)/2, -n), (sympy.S.Half,), -b*x/a)/(sqrt(x)*(1 + b*x/a)**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_748():
    f = (a + b*x)**n/x**(sympy.S(5)/2)
    F = -2*(a + b*x)**n*hyper((sympy.S(-3)/2, -n), (sympy.S(-1)/2,), -b*x/a)/(3*x**(sympy.S(3)/2)*(1 + b*x/a)**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_749():
    f = (b*x)**m*(d*x + 2)**n
    F = 2**n*(b*x)**(m + 1)*hyper((-n, m + 1), (m + 2,), -d*x/2)/(b*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_750():
    f = (b*x)**m*(-b*c*x + c)**n
    F = -(-b*c*x + c)**(n + 1)*hyper((-m, n + 1), (n + 2,), -b*x + 1)/(b*c*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_751():
    f = (b*x)**m*(c + d*x)**n
    F = (b*x)**(m + 1)*(c + d*x)**n*hyper((-n, m + 1), (m + 2,), -d*x/c)/(b*(1 + d*x/c)**n*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_752():
    f = x**(n - 1)*(a + b*x)**(-n - 1)
    F = x**n/(a*n*(a + b*x)**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_753():
    f = x**(-n - 3)*(a + b*x)**n
    F = -x**(-n - 2)*(a + b*x)**(n + 1)/(a*(n + 2)) + b*x**(-n - 1)*(a + b*x)**(n + 1)/(a**2*(n + 1)*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_754():
    f = x**(-n - 3)*(a + b*x)**n
    F = -x**(-n - 2)*(a + b*x)**(n + 1)/(a*(n + 2)) + b*x**(-n - 1)*(a + b*x)**(n + 1)/(a**2*(n + 1)*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_755():
    f = x**3*sqrt(c*x**2)*(a + b*x)
    F = a*x**4*sqrt(c*x**2)/5 + b*x**5*sqrt(c*x**2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_756():
    f = x**2*sqrt(c*x**2)*(a + b*x)
    F = a*x**3*sqrt(c*x**2)/4 + b*x**4*sqrt(c*x**2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_757():
    f = x*sqrt(c*x**2)*(a + b*x)
    F = a*x**2*sqrt(c*x**2)/3 + b*x**3*sqrt(c*x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_758():
    f = sqrt(c*x**2)*(a + b*x)
    F = a*x*sqrt(c*x**2)/2 + b*x**2*sqrt(c*x**2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_759():
    f = sqrt(c*x**2)*(a + b*x)/x
    F = a*sqrt(c*x**2) + b*x*sqrt(c*x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_760():
    f = sqrt(c*x**2)*(a + b*x)/x**2
    F = a*sqrt(c*x**2)*log(x)/x + b*sqrt(c*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_761():
    f = sqrt(c*x**2)*(a + b*x)/x**3
    F = -a*sqrt(c*x**2)/x**2 + b*sqrt(c*x**2)*log(x)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_762():
    f = sqrt(c*x**2)*(a + b*x)/x**4
    F = -sqrt(c*x**2)*(a + b*x)**2/(2*a*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_763():
    f = x**3*(c*x**2)**(sympy.S(3)/2)*(a + b*x)
    F = a*c*x**6*sqrt(c*x**2)/7 + b*c*x**7*sqrt(c*x**2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_764():
    f = x**2*(c*x**2)**(sympy.S(3)/2)*(a + b*x)
    F = a*c*x**5*sqrt(c*x**2)/6 + b*c*x**6*sqrt(c*x**2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_765():
    f = x*(c*x**2)**(sympy.S(3)/2)*(a + b*x)
    F = a*c*x**4*sqrt(c*x**2)/5 + b*c*x**5*sqrt(c*x**2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_766():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)
    F = a*c*x**3*sqrt(c*x**2)/4 + b*c*x**4*sqrt(c*x**2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_767():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)/x
    F = a*c*x**2*sqrt(c*x**2)/3 + b*c*x**3*sqrt(c*x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_768():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)/x**2
    F = a*c*x*sqrt(c*x**2)/2 + b*c*x**2*sqrt(c*x**2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_769():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)/x**3
    F = a*c*sqrt(c*x**2) + b*c*x*sqrt(c*x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_770():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)/x**4
    F = a*c*sqrt(c*x**2)*log(x)/x + b*c*sqrt(c*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_771():
    f = x**3*(c*x**2)**(sympy.S(5)/2)*(a + b*x)
    F = a*c**2*x**8*sqrt(c*x**2)/9 + b*c**2*x**9*sqrt(c*x**2)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_772():
    f = x**2*(c*x**2)**(sympy.S(5)/2)*(a + b*x)
    F = a*c**2*x**7*sqrt(c*x**2)/8 + b*c**2*x**8*sqrt(c*x**2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_773():
    f = x*(c*x**2)**(sympy.S(5)/2)*(a + b*x)
    F = a*c**2*x**6*sqrt(c*x**2)/7 + b*c**2*x**7*sqrt(c*x**2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_774():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)
    F = a*c**2*x**5*sqrt(c*x**2)/6 + b*c**2*x**6*sqrt(c*x**2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_775():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)/x
    F = a*c**2*x**4*sqrt(c*x**2)/5 + b*c**2*x**5*sqrt(c*x**2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_776():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)/x**2
    F = a*c**2*x**3*sqrt(c*x**2)/4 + b*c**2*x**4*sqrt(c*x**2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_777():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)/x**3
    F = a*c**2*x**2*sqrt(c*x**2)/3 + b*c**2*x**3*sqrt(c*x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_778():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)/x**4
    F = a*c**2*x*sqrt(c*x**2)/2 + b*c**2*x**2*sqrt(c*x**2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_779():
    f = x**3*(a + b*x)/sqrt(c*x**2)
    F = a*x**4/(3*sqrt(c*x**2)) + b*x**5/(4*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_780():
    f = x**2*(a + b*x)/sqrt(c*x**2)
    F = a*x**3/(2*sqrt(c*x**2)) + b*x**4/(3*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_781():
    f = x*(a + b*x)/sqrt(c*x**2)
    F = a*x**2/sqrt(c*x**2) + b*x**3/(2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_782():
    f = (a + b*x)/sqrt(c*x**2)
    F = a*x*log(x)/sqrt(c*x**2) + b*x**2/sqrt(c*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_783():
    f = (a + b*x)/(x*sqrt(c*x**2))
    F = -a/sqrt(c*x**2) + b*x*log(x)/sqrt(c*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_784():
    f = (a + b*x)/(x**2*sqrt(c*x**2))
    F = -(a + b*x)**2/(2*a*x*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_785():
    f = (a + b*x)/(x**3*sqrt(c*x**2))
    F = -a/(3*x**2*sqrt(c*x**2)) - b/(2*x*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_786():
    f = (a + b*x)/(x**4*sqrt(c*x**2))
    F = -a/(4*x**3*sqrt(c*x**2)) - b/(3*x**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_787():
    f = x**3*(a + b*x)/(c*x**2)**(sympy.S(3)/2)
    F = a*x**2/(c*sqrt(c*x**2)) + b*x**3/(2*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_788():
    f = x**2*(a + b*x)/(c*x**2)**(sympy.S(3)/2)
    F = a*x*log(x)/(c*sqrt(c*x**2)) + b*x**2/(c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_789():
    f = x*(a + b*x)/(c*x**2)**(sympy.S(3)/2)
    F = -a/(c*sqrt(c*x**2)) + b*x*log(x)/(c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_790():
    f = (a + b*x)/(c*x**2)**(sympy.S(3)/2)
    F = -(a + b*x)**2/(2*a*c*x*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_791():
    f = (a + b*x)/(x*(c*x**2)**(sympy.S(3)/2))
    F = -a/(3*c*x**2*sqrt(c*x**2)) - b/(2*c*x*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_792():
    f = (a + b*x)/(x**2*(c*x**2)**(sympy.S(3)/2))
    F = -a/(4*c*x**3*sqrt(c*x**2)) - b/(3*c*x**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_793():
    f = (a + b*x)/(x**3*(c*x**2)**(sympy.S(3)/2))
    F = -a/(5*c*x**4*sqrt(c*x**2)) - b/(4*c*x**3*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_794():
    f = (a + b*x)/(x**4*(c*x**2)**(sympy.S(3)/2))
    F = -a/(6*c*x**5*sqrt(c*x**2)) - b/(5*c*x**4*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_795():
    f = x**3*(a + b*x)/(c*x**2)**(sympy.S(5)/2)
    F = -a/(c**2*sqrt(c*x**2)) + b*x*log(x)/(c**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_796():
    f = x**2*(a + b*x)/(c*x**2)**(sympy.S(5)/2)
    F = -(a + b*x)**2/(2*a*c**2*x*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_797():
    f = x*(a + b*x)/(c*x**2)**(sympy.S(5)/2)
    F = -a/(3*c**2*x**2*sqrt(c*x**2)) - b/(2*c**2*x*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_798():
    f = (a + b*x)/(c*x**2)**(sympy.S(5)/2)
    F = -a/(4*c**2*x**3*sqrt(c*x**2)) - b/(3*c**2*x**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_799():
    f = (a + b*x)/(x*(c*x**2)**(sympy.S(5)/2))
    F = -a/(5*c**2*x**4*sqrt(c*x**2)) - b/(4*c**2*x**3*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_800():
    f = (a + b*x)/(x**2*(c*x**2)**(sympy.S(5)/2))
    F = -a/(6*c**2*x**5*sqrt(c*x**2)) - b/(5*c**2*x**4*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_801():
    f = (a + b*x)/(x**3*(c*x**2)**(sympy.S(5)/2))
    F = -a/(7*c**2*x**6*sqrt(c*x**2)) - b/(6*c**2*x**5*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_802():
    f = (a + b*x)/(x**4*(c*x**2)**(sympy.S(5)/2))
    F = -a/(8*c**2*x**7*sqrt(c*x**2)) - b/(7*c**2*x**6*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_803():
    f = x**3*sqrt(c*x**2)*(a + b*x)**2
    F = a**2*x**4*sqrt(c*x**2)/5 + a*b*x**5*sqrt(c*x**2)/3 + b**2*x**6*sqrt(c*x**2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_804():
    f = x**2*sqrt(c*x**2)*(a + b*x)**2
    F = a**2*x**3*sqrt(c*x**2)/4 + 2*a*b*x**4*sqrt(c*x**2)/5 + b**2*x**5*sqrt(c*x**2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_805():
    f = x*sqrt(c*x**2)*(a + b*x)**2
    F = a**2*x**2*sqrt(c*x**2)/3 + a*b*x**3*sqrt(c*x**2)/2 + b**2*x**4*sqrt(c*x**2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_806():
    f = sqrt(c*x**2)*(a + b*x)**2
    F = a**2*x*sqrt(c*x**2)/2 + 2*a*b*x**2*sqrt(c*x**2)/3 + b**2*x**3*sqrt(c*x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_807():
    f = sqrt(c*x**2)*(a + b*x)**2/x
    F = sqrt(c*x**2)*(a + b*x)**3/(3*b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_808():
    f = sqrt(c*x**2)*(a + b*x)**2/x**2
    F = a**2*sqrt(c*x**2)*log(x)/x + 2*a*b*sqrt(c*x**2) + b**2*x*sqrt(c*x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_809():
    f = sqrt(c*x**2)*(a + b*x)**2/x**3
    F = -a**2*sqrt(c*x**2)/x**2 + 2*a*b*sqrt(c*x**2)*log(x)/x + b**2*sqrt(c*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_810():
    f = sqrt(c*x**2)*(a + b*x)**2/x**4
    F = -a**2*sqrt(c*x**2)/(2*x**3) - 2*a*b*sqrt(c*x**2)/x**2 + b**2*sqrt(c*x**2)*log(x)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_811():
    f = x**3*(c*x**2)**(sympy.S(3)/2)*(a + b*x)**2
    F = a**2*c*x**6*sqrt(c*x**2)/7 + a*b*c*x**7*sqrt(c*x**2)/4 + b**2*c*x**8*sqrt(c*x**2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_812():
    f = x**2*(c*x**2)**(sympy.S(3)/2)*(a + b*x)**2
    F = a**2*c*x**5*sqrt(c*x**2)/6 + 2*a*b*c*x**6*sqrt(c*x**2)/7 + b**2*c*x**7*sqrt(c*x**2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_813():
    f = x*(c*x**2)**(sympy.S(3)/2)*(a + b*x)**2
    F = a**2*c*x**4*sqrt(c*x**2)/5 + a*b*c*x**5*sqrt(c*x**2)/3 + b**2*c*x**6*sqrt(c*x**2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_814():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)**2
    F = a**2*c*x**3*sqrt(c*x**2)/4 + 2*a*b*c*x**4*sqrt(c*x**2)/5 + b**2*c*x**5*sqrt(c*x**2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_815():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)**2/x
    F = a**2*c*x**2*sqrt(c*x**2)/3 + a*b*c*x**3*sqrt(c*x**2)/2 + b**2*c*x**4*sqrt(c*x**2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_816():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)**2/x**2
    F = a**2*c*x*sqrt(c*x**2)/2 + 2*a*b*c*x**2*sqrt(c*x**2)/3 + b**2*c*x**3*sqrt(c*x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_817():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)**2/x**3
    F = c*sqrt(c*x**2)*(a + b*x)**3/(3*b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_818():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)**2/x**4
    F = a**2*c*sqrt(c*x**2)*log(x)/x + 2*a*b*c*sqrt(c*x**2) + b**2*c*x*sqrt(c*x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_819():
    f = x*(c*x**2)**(sympy.S(5)/2)*(a + b*x)**2
    F = a**2*c**2*x**6*sqrt(c*x**2)/7 + a*b*c**2*x**7*sqrt(c*x**2)/4 + b**2*c**2*x**8*sqrt(c*x**2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_820():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)**2
    F = a**2*c**2*x**5*sqrt(c*x**2)/6 + 2*a*b*c**2*x**6*sqrt(c*x**2)/7 + b**2*c**2*x**7*sqrt(c*x**2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_821():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)**2/x
    F = a**2*c**2*x**4*sqrt(c*x**2)/5 + a*b*c**2*x**5*sqrt(c*x**2)/3 + b**2*c**2*x**6*sqrt(c*x**2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_822():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)**2/x**2
    F = a**2*c**2*x**3*sqrt(c*x**2)/4 + 2*a*b*c**2*x**4*sqrt(c*x**2)/5 + b**2*c**2*x**5*sqrt(c*x**2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_823():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)**2/x**3
    F = a**2*c**2*x**2*sqrt(c*x**2)/3 + a*b*c**2*x**3*sqrt(c*x**2)/2 + b**2*c**2*x**4*sqrt(c*x**2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_824():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)**2/x**4
    F = a**2*c**2*x*sqrt(c*x**2)/2 + 2*a*b*c**2*x**2*sqrt(c*x**2)/3 + b**2*c**2*x**3*sqrt(c*x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_825():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)**2/x**5
    F = c**2*sqrt(c*x**2)*(a + b*x)**3/(3*b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_826():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)**2/x**6
    F = a**2*c**2*sqrt(c*x**2)*log(x)/x + 2*a*b*c**2*sqrt(c*x**2) + b**2*c**2*x*sqrt(c*x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_827():
    f = x**3*(a + b*x)**2/sqrt(c*x**2)
    F = a**2*x**4/(3*sqrt(c*x**2)) + a*b*x**5/(2*sqrt(c*x**2)) + b**2*x**6/(5*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_828():
    f = x**2*(a + b*x)**2/sqrt(c*x**2)
    F = a**2*x**3/(2*sqrt(c*x**2)) + 2*a*b*x**4/(3*sqrt(c*x**2)) + b**2*x**5/(4*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_829():
    f = x*(a + b*x)**2/sqrt(c*x**2)
    F = x*(a + b*x)**3/(3*b*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_830():
    f = (a + b*x)**2/sqrt(c*x**2)
    F = a**2*x*log(x)/sqrt(c*x**2) + 2*a*b*x**2/sqrt(c*x**2) + b**2*x**3/(2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_831():
    f = (a + b*x)**2/(x*sqrt(c*x**2))
    F = -a**2/sqrt(c*x**2) + 2*a*b*x*log(x)/sqrt(c*x**2) + b**2*x**2/sqrt(c*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_832():
    f = (a + b*x)**2/(x**2*sqrt(c*x**2))
    F = -a**2/(2*x*sqrt(c*x**2)) - 2*a*b/sqrt(c*x**2) + b**2*x*log(x)/sqrt(c*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_833():
    f = (a + b*x)**2/(x**3*sqrt(c*x**2))
    F = -(a + b*x)**3/(3*a*x**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_834():
    f = (a + b*x)**2/(x**4*sqrt(c*x**2))
    F = -a**2/(4*x**3*sqrt(c*x**2)) - 2*a*b/(3*x**2*sqrt(c*x**2)) - b**2/(2*x*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_835():
    f = x**3*(a + b*x)**2/(c*x**2)**(sympy.S(3)/2)
    F = x*(a + b*x)**3/(3*b*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_836():
    f = x**2*(a + b*x)**2/(c*x**2)**(sympy.S(3)/2)
    F = a**2*x*log(x)/(c*sqrt(c*x**2)) + 2*a*b*x**2/(c*sqrt(c*x**2)) + b**2*x**3/(2*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_837():
    f = x*(a + b*x)**2/(c*x**2)**(sympy.S(3)/2)
    F = -a**2/(c*sqrt(c*x**2)) + 2*a*b*x*log(x)/(c*sqrt(c*x**2)) + b**2*x**2/(c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_838():
    f = (a + b*x)**2/(c*x**2)**(sympy.S(3)/2)
    F = -a**2/(2*c*x*sqrt(c*x**2)) - 2*a*b/(c*sqrt(c*x**2)) + b**2*x*log(x)/(c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_839():
    f = (a + b*x)**2/(x*(c*x**2)**(sympy.S(3)/2))
    F = -(a + b*x)**3/(3*a*c*x**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_840():
    f = (a + b*x)**2/(x**2*(c*x**2)**(sympy.S(3)/2))
    F = -a**2/(4*c*x**3*sqrt(c*x**2)) - 2*a*b/(3*c*x**2*sqrt(c*x**2)) - b**2/(2*c*x*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_841():
    f = (a + b*x)**2/(x**3*(c*x**2)**(sympy.S(3)/2))
    F = -a**2/(5*c*x**4*sqrt(c*x**2)) - a*b/(2*c*x**3*sqrt(c*x**2)) - b**2/(3*c*x**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_842():
    f = (a + b*x)**2/(x**4*(c*x**2)**(sympy.S(3)/2))
    F = -a**2/(6*c*x**5*sqrt(c*x**2)) - 2*a*b/(5*c*x**4*sqrt(c*x**2)) - b**2/(4*c*x**3*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_843():
    f = x**3*(a + b*x)**2/(c*x**2)**(sympy.S(5)/2)
    F = -a**2/(c**2*sqrt(c*x**2)) + 2*a*b*x*log(x)/(c**2*sqrt(c*x**2)) + b**2*x**2/(c**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_844():
    f = x**2*(a + b*x)**2/(c*x**2)**(sympy.S(5)/2)
    F = -a**2/(2*c**2*x*sqrt(c*x**2)) - 2*a*b/(c**2*sqrt(c*x**2)) + b**2*x*log(x)/(c**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_845():
    f = x*(a + b*x)**2/(c*x**2)**(sympy.S(5)/2)
    F = -(a + b*x)**3/(3*a*c**2*x**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_846():
    f = (a + b*x)**2/(c*x**2)**(sympy.S(5)/2)
    F = -a**2/(4*c**2*x**3*sqrt(c*x**2)) - 2*a*b/(3*c**2*x**2*sqrt(c*x**2)) - b**2/(2*c**2*x*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_847():
    f = (a + b*x)**2/(x*(c*x**2)**(sympy.S(5)/2))
    F = -a**2/(5*c**2*x**4*sqrt(c*x**2)) - a*b/(2*c**2*x**3*sqrt(c*x**2)) - b**2/(3*c**2*x**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_848():
    f = (a + b*x)**2/(x**2*(c*x**2)**(sympy.S(5)/2))
    F = -a**2/(6*c**2*x**5*sqrt(c*x**2)) - 2*a*b/(5*c**2*x**4*sqrt(c*x**2)) - b**2/(4*c**2*x**3*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_849():
    f = (a + b*x)**2/(x**3*(c*x**2)**(sympy.S(5)/2))
    F = -a**2/(7*c**2*x**6*sqrt(c*x**2)) - a*b/(3*c**2*x**5*sqrt(c*x**2)) - b**2/(5*c**2*x**4*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_850():
    f = (a + b*x)**2/(x**4*(c*x**2)**(sympy.S(5)/2))
    F = -a**2/(8*c**2*x**7*sqrt(c*x**2)) - 2*a*b/(7*c**2*x**6*sqrt(c*x**2)) - b**2/(6*c**2*x**5*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_851():
    f = x**3*sqrt(c*x**2)/(a + b*x)
    F = a**4*sqrt(c*x**2)*log(a + b*x)/(b**5*x) - a**3*sqrt(c*x**2)/b**4 + a**2*x*sqrt(c*x**2)/(2*b**3) - a*x**2*sqrt(c*x**2)/(3*b**2) + x**3*sqrt(c*x**2)/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_852():
    f = x**2*sqrt(c*x**2)/(a + b*x)
    F = -a**3*sqrt(c*x**2)*log(a + b*x)/(b**4*x) + a**2*sqrt(c*x**2)/b**3 - a*x*sqrt(c*x**2)/(2*b**2) + x**2*sqrt(c*x**2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_853():
    f = x*sqrt(c*x**2)/(a + b*x)
    F = a**2*sqrt(c*x**2)*log(a + b*x)/(b**3*x) - a*sqrt(c*x**2)/b**2 + x*sqrt(c*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_854():
    f = sqrt(c*x**2)/(a + b*x)
    F = -a*sqrt(c*x**2)*log(a + b*x)/(b**2*x) + sqrt(c*x**2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_855():
    f = sqrt(c*x**2)/(x*(a + b*x))
    F = sqrt(c*x**2)*log(a + b*x)/(b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_856():
    f = sqrt(c*x**2)/(x**2*(a + b*x))
    F = sqrt(c*x**2)*log(x)/(a*x) - sqrt(c*x**2)*log(a + b*x)/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_857():
    f = sqrt(c*x**2)/(x**3*(a + b*x))
    F = -sqrt(c*x**2)/(a*x**2) - b*sqrt(c*x**2)*log(x)/(a**2*x) + b*sqrt(c*x**2)*log(a + b*x)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_858():
    f = sqrt(c*x**2)/(x**4*(a + b*x))
    F = -sqrt(c*x**2)/(2*a*x**3) + b*sqrt(c*x**2)/(a**2*x**2) + b**2*sqrt(c*x**2)*log(x)/(a**3*x) - b**2*sqrt(c*x**2)*log(a + b*x)/(a**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_859():
    f = x*(c*x**2)**(sympy.S(3)/2)/(a + b*x)
    F = a**4*c*sqrt(c*x**2)*log(a + b*x)/(b**5*x) - a**3*c*sqrt(c*x**2)/b**4 + a**2*c*x*sqrt(c*x**2)/(2*b**3) - a*c*x**2*sqrt(c*x**2)/(3*b**2) + c*x**3*sqrt(c*x**2)/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_860():
    f = (c*x**2)**(sympy.S(3)/2)/(a + b*x)
    F = -a**3*c*sqrt(c*x**2)*log(a + b*x)/(b**4*x) + a**2*c*sqrt(c*x**2)/b**3 - a*c*x*sqrt(c*x**2)/(2*b**2) + c*x**2*sqrt(c*x**2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_861():
    f = (c*x**2)**(sympy.S(3)/2)/(x*(a + b*x))
    F = a**2*c*sqrt(c*x**2)*log(a + b*x)/(b**3*x) - a*c*sqrt(c*x**2)/b**2 + c*x*sqrt(c*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_862():
    f = (c*x**2)**(sympy.S(3)/2)/(x**2*(a + b*x))
    F = -a*c*sqrt(c*x**2)*log(a + b*x)/(b**2*x) + c*sqrt(c*x**2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_863():
    f = (c*x**2)**(sympy.S(3)/2)/(x**3*(a + b*x))
    F = c*sqrt(c*x**2)*log(a + b*x)/(b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_864():
    f = (c*x**2)**(sympy.S(3)/2)/(x**4*(a + b*x))
    F = c*sqrt(c*x**2)*log(x)/(a*x) - c*sqrt(c*x**2)*log(a + b*x)/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_865():
    f = (c*x**2)**(sympy.S(3)/2)/(x**5*(a + b*x))
    F = -c*sqrt(c*x**2)/(a*x**2) - b*c*sqrt(c*x**2)*log(x)/(a**2*x) + b*c*sqrt(c*x**2)*log(a + b*x)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_866():
    f = (c*x**2)**(sympy.S(3)/2)/(x**6*(a + b*x))
    F = -c*sqrt(c*x**2)/(2*a*x**3) + b*c*sqrt(c*x**2)/(a**2*x**2) + b**2*c*sqrt(c*x**2)*log(x)/(a**3*x) - b**2*c*sqrt(c*x**2)*log(a + b*x)/(a**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_867():
    f = (c*x**2)**(sympy.S(3)/2)/(x**7*(a + b*x))
    F = -c*sqrt(c*x**2)/(3*a*x**4) + b*c*sqrt(c*x**2)/(2*a**2*x**3) - b**2*c*sqrt(c*x**2)/(a**3*x**2) - b**3*c*sqrt(c*x**2)*log(x)/(a**4*x) + b**3*c*sqrt(c*x**2)*log(a + b*x)/(a**4*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_868():
    f = (c*x**2)**(sympy.S(5)/2)/(a + b*x)
    F = -a**5*c**2*sqrt(c*x**2)*log(a + b*x)/(b**6*x) + a**4*c**2*sqrt(c*x**2)/b**5 - a**3*c**2*x*sqrt(c*x**2)/(2*b**4) + a**2*c**2*x**2*sqrt(c*x**2)/(3*b**3) - a*c**2*x**3*sqrt(c*x**2)/(4*b**2) + c**2*x**4*sqrt(c*x**2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_869():
    f = (c*x**2)**(sympy.S(5)/2)/(x*(a + b*x))
    F = a**4*c**2*sqrt(c*x**2)*log(a + b*x)/(b**5*x) - a**3*c**2*sqrt(c*x**2)/b**4 + a**2*c**2*x*sqrt(c*x**2)/(2*b**3) - a*c**2*x**2*sqrt(c*x**2)/(3*b**2) + c**2*x**3*sqrt(c*x**2)/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_870():
    f = (c*x**2)**(sympy.S(5)/2)/(x**2*(a + b*x))
    F = -a**3*c**2*sqrt(c*x**2)*log(a + b*x)/(b**4*x) + a**2*c**2*sqrt(c*x**2)/b**3 - a*c**2*x*sqrt(c*x**2)/(2*b**2) + c**2*x**2*sqrt(c*x**2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_871():
    f = (c*x**2)**(sympy.S(5)/2)/(x**3*(a + b*x))
    F = a**2*c**2*sqrt(c*x**2)*log(a + b*x)/(b**3*x) - a*c**2*sqrt(c*x**2)/b**2 + c**2*x*sqrt(c*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_872():
    f = (c*x**2)**(sympy.S(5)/2)/(x**4*(a + b*x))
    F = -a*c**2*sqrt(c*x**2)*log(a + b*x)/(b**2*x) + c**2*sqrt(c*x**2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_873():
    f = (c*x**2)**(sympy.S(5)/2)/(x**5*(a + b*x))
    F = c**2*sqrt(c*x**2)*log(a + b*x)/(b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_874():
    f = (c*x**2)**(sympy.S(5)/2)/(x**6*(a + b*x))
    F = c**2*sqrt(c*x**2)*log(x)/(a*x) - c**2*sqrt(c*x**2)*log(a + b*x)/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_875():
    f = (c*x**2)**(sympy.S(5)/2)/(x**7*(a + b*x))
    F = -c**2*sqrt(c*x**2)/(a*x**2) - b*c**2*sqrt(c*x**2)*log(x)/(a**2*x) + b*c**2*sqrt(c*x**2)*log(a + b*x)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_876():
    f = x**4/(sqrt(c*x**2)*(a + b*x))
    F = -a**3*x*log(a + b*x)/(b**4*sqrt(c*x**2)) + a**2*x**2/(b**3*sqrt(c*x**2)) - a*x**3/(2*b**2*sqrt(c*x**2)) + x**4/(3*b*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_877():
    f = x**3/(sqrt(c*x**2)*(a + b*x))
    F = a**2*x*log(a + b*x)/(b**3*sqrt(c*x**2)) - a*x**2/(b**2*sqrt(c*x**2)) + x**3/(2*b*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_878():
    f = x**2/(sqrt(c*x**2)*(a + b*x))
    F = -a*x*log(a + b*x)/(b**2*sqrt(c*x**2)) + x**2/(b*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_879():
    f = x/(sqrt(c*x**2)*(a + b*x))
    F = x*log(a + b*x)/(b*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_880():
    f = 1/(sqrt(c*x**2)*(a + b*x))
    F = x*log(x)/(a*sqrt(c*x**2)) - x*log(a + b*x)/(a*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_881():
    f = 1/(x*sqrt(c*x**2)*(a + b*x))
    F = -1/(a*sqrt(c*x**2)) - b*x*log(x)/(a**2*sqrt(c*x**2)) + b*x*log(a + b*x)/(a**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_882():
    f = 1/(x**2*sqrt(c*x**2)*(a + b*x))
    F = -1/(2*a*x*sqrt(c*x**2)) + b/(a**2*sqrt(c*x**2)) + b**2*x*log(x)/(a**3*sqrt(c*x**2)) - b**2*x*log(a + b*x)/(a**3*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_883():
    f = 1/(x**3*sqrt(c*x**2)*(a + b*x))
    F = -1/(3*a*x**2*sqrt(c*x**2)) + b/(2*a**2*x*sqrt(c*x**2)) - b**2/(a**3*sqrt(c*x**2)) - b**3*x*log(x)/(a**4*sqrt(c*x**2)) + b**3*x*log(a + b*x)/(a**4*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_884():
    f = x**6/((c*x**2)**(sympy.S(3)/2)*(a + b*x))
    F = -a**3*x*log(a + b*x)/(b**4*c*sqrt(c*x**2)) + a**2*x**2/(b**3*c*sqrt(c*x**2)) - a*x**3/(2*b**2*c*sqrt(c*x**2)) + x**4/(3*b*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_885():
    f = x**5/((c*x**2)**(sympy.S(3)/2)*(a + b*x))
    F = a**2*x*log(a + b*x)/(b**3*c*sqrt(c*x**2)) - a*x**2/(b**2*c*sqrt(c*x**2)) + x**3/(2*b*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_886():
    f = x**4/((c*x**2)**(sympy.S(3)/2)*(a + b*x))
    F = -a*x*log(a + b*x)/(b**2*c*sqrt(c*x**2)) + x**2/(b*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_887():
    f = x**3/((c*x**2)**(sympy.S(3)/2)*(a + b*x))
    F = x*log(a + b*x)/(b*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_888():
    f = x**2/((c*x**2)**(sympy.S(3)/2)*(a + b*x))
    F = x*log(x)/(a*c*sqrt(c*x**2)) - x*log(a + b*x)/(a*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_889():
    f = x/((c*x**2)**(sympy.S(3)/2)*(a + b*x))
    F = -1/(a*c*sqrt(c*x**2)) - b*x*log(x)/(a**2*c*sqrt(c*x**2)) + b*x*log(a + b*x)/(a**2*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_890():
    f = 1/((c*x**2)**(sympy.S(3)/2)*(a + b*x))
    F = -1/(2*a*c*x*sqrt(c*x**2)) + b/(a**2*c*sqrt(c*x**2)) + b**2*x*log(x)/(a**3*c*sqrt(c*x**2)) - b**2*x*log(a + b*x)/(a**3*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_891():
    f = 1/(x*(c*x**2)**(sympy.S(3)/2)*(a + b*x))
    F = -1/(3*a*c*x**2*sqrt(c*x**2)) + b/(2*a**2*c*x*sqrt(c*x**2)) - b**2/(a**3*c*sqrt(c*x**2)) - b**3*x*log(x)/(a**4*c*sqrt(c*x**2)) + b**3*x*log(a + b*x)/(a**4*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_892():
    f = x**3*sqrt(c*x**2)/(a + b*x)**2
    F = -a**4*sqrt(c*x**2)/(b**5*x*(a + b*x)) - 4*a**3*sqrt(c*x**2)*log(a + b*x)/(b**5*x) + 3*a**2*sqrt(c*x**2)/b**4 - a*x*sqrt(c*x**2)/b**3 + x**2*sqrt(c*x**2)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_893():
    f = x**2*sqrt(c*x**2)/(a + b*x)**2
    F = a**3*sqrt(c*x**2)/(b**4*x*(a + b*x)) + 3*a**2*sqrt(c*x**2)*log(a + b*x)/(b**4*x) - 2*a*sqrt(c*x**2)/b**3 + x*sqrt(c*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_894():
    f = x*sqrt(c*x**2)/(a + b*x)**2
    F = -a**2*sqrt(c*x**2)/(b**3*x*(a + b*x)) - 2*a*sqrt(c*x**2)*log(a + b*x)/(b**3*x) + sqrt(c*x**2)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_895():
    f = sqrt(c*x**2)/(a + b*x)**2
    F = a*sqrt(c*x**2)/(b**2*x*(a + b*x)) + sqrt(c*x**2)*log(a + b*x)/(b**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_896():
    f = sqrt(c*x**2)/(x*(a + b*x)**2)
    F = -sqrt(c*x**2)/(b*x*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_897():
    f = sqrt(c*x**2)/(x**2*(a + b*x)**2)
    F = sqrt(c*x**2)/(a*x*(a + b*x)) + sqrt(c*x**2)*log(x)/(a**2*x) - sqrt(c*x**2)*log(a + b*x)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_898():
    f = sqrt(c*x**2)/(x**3*(a + b*x)**2)
    F = -b*sqrt(c*x**2)/(a**2*x*(a + b*x)) - sqrt(c*x**2)/(a**2*x**2) - 2*b*sqrt(c*x**2)*log(x)/(a**3*x) + 2*b*sqrt(c*x**2)*log(a + b*x)/(a**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_899():
    f = sqrt(c*x**2)/(x**4*(a + b*x)**2)
    F = -sqrt(c*x**2)/(2*a**2*x**3) + b**2*sqrt(c*x**2)/(a**3*x*(a + b*x)) + 2*b*sqrt(c*x**2)/(a**3*x**2) + 3*b**2*sqrt(c*x**2)*log(x)/(a**4*x) - 3*b**2*sqrt(c*x**2)*log(a + b*x)/(a**4*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_900():
    f = x*(c*x**2)**(sympy.S(3)/2)/(a + b*x)**2
    F = -a**4*c*sqrt(c*x**2)/(b**5*x*(a + b*x)) - 4*a**3*c*sqrt(c*x**2)*log(a + b*x)/(b**5*x) + 3*a**2*c*sqrt(c*x**2)/b**4 - a*c*x*sqrt(c*x**2)/b**3 + c*x**2*sqrt(c*x**2)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_901():
    f = (c*x**2)**(sympy.S(3)/2)/(a + b*x)**2
    F = a**3*c*sqrt(c*x**2)/(b**4*x*(a + b*x)) + 3*a**2*c*sqrt(c*x**2)*log(a + b*x)/(b**4*x) - 2*a*c*sqrt(c*x**2)/b**3 + c*x*sqrt(c*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_902():
    f = (c*x**2)**(sympy.S(3)/2)/(x*(a + b*x)**2)
    F = -a**2*c*sqrt(c*x**2)/(b**3*x*(a + b*x)) - 2*a*c*sqrt(c*x**2)*log(a + b*x)/(b**3*x) + c*sqrt(c*x**2)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_903():
    f = (c*x**2)**(sympy.S(3)/2)/(x**2*(a + b*x)**2)
    F = a*c*sqrt(c*x**2)/(b**2*x*(a + b*x)) + c*sqrt(c*x**2)*log(a + b*x)/(b**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_904():
    f = (c*x**2)**(sympy.S(3)/2)/(x**3*(a + b*x)**2)
    F = -c*sqrt(c*x**2)/(b*x*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_905():
    f = (c*x**2)**(sympy.S(3)/2)/(x**4*(a + b*x)**2)
    F = c*sqrt(c*x**2)/(a*x*(a + b*x)) + c*sqrt(c*x**2)*log(x)/(a**2*x) - c*sqrt(c*x**2)*log(a + b*x)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_906():
    f = (c*x**2)**(sympy.S(3)/2)/(x**5*(a + b*x)**2)
    F = -b*c*sqrt(c*x**2)/(a**2*x*(a + b*x)) - c*sqrt(c*x**2)/(a**2*x**2) - 2*b*c*sqrt(c*x**2)*log(x)/(a**3*x) + 2*b*c*sqrt(c*x**2)*log(a + b*x)/(a**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_907():
    f = (c*x**2)**(sympy.S(3)/2)/(x**6*(a + b*x)**2)
    F = -c*sqrt(c*x**2)/(2*a**2*x**3) + b**2*c*sqrt(c*x**2)/(a**3*x*(a + b*x)) + 2*b*c*sqrt(c*x**2)/(a**3*x**2) + 3*b**2*c*sqrt(c*x**2)*log(x)/(a**4*x) - 3*b**2*c*sqrt(c*x**2)*log(a + b*x)/(a**4*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_908():
    f = x**5/(sqrt(c*x**2)*(a + b*x)**2)
    F = -a**4*x/(b**5*sqrt(c*x**2)*(a + b*x)) - 4*a**3*x*log(a + b*x)/(b**5*sqrt(c*x**2)) + 3*a**2*x**2/(b**4*sqrt(c*x**2)) - a*x**3/(b**3*sqrt(c*x**2)) + x**4/(3*b**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_909():
    f = x**4/(sqrt(c*x**2)*(a + b*x)**2)
    F = a**3*x/(b**4*sqrt(c*x**2)*(a + b*x)) + 3*a**2*x*log(a + b*x)/(b**4*sqrt(c*x**2)) - 2*a*x**2/(b**3*sqrt(c*x**2)) + x**3/(2*b**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_910():
    f = x**3/(sqrt(c*x**2)*(a + b*x)**2)
    F = -a**2*x/(b**3*sqrt(c*x**2)*(a + b*x)) - 2*a*x*log(a + b*x)/(b**3*sqrt(c*x**2)) + x**2/(b**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_911():
    f = x**2/(sqrt(c*x**2)*(a + b*x)**2)
    F = a*x/(b**2*sqrt(c*x**2)*(a + b*x)) + x*log(a + b*x)/(b**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_912():
    f = x/(sqrt(c*x**2)*(a + b*x)**2)
    F = -x/(b*sqrt(c*x**2)*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_913():
    f = 1/(sqrt(c*x**2)*(a + b*x)**2)
    F = x/(a*sqrt(c*x**2)*(a + b*x)) + x*log(x)/(a**2*sqrt(c*x**2)) - x*log(a + b*x)/(a**2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_914():
    f = 1/(x*sqrt(c*x**2)*(a + b*x)**2)
    F = -b*x/(a**2*sqrt(c*x**2)*(a + b*x)) - 1/(a**2*sqrt(c*x**2)) - 2*b*x*log(x)/(a**3*sqrt(c*x**2)) + 2*b*x*log(a + b*x)/(a**3*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_915():
    f = 1/(x**2*sqrt(c*x**2)*(a + b*x)**2)
    F = -1/(2*a**2*x*sqrt(c*x**2)) + b**2*x/(a**3*sqrt(c*x**2)*(a + b*x)) + 2*b/(a**3*sqrt(c*x**2)) + 3*b**2*x*log(x)/(a**4*sqrt(c*x**2)) - 3*b**2*x*log(a + b*x)/(a**4*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_916():
    f = x**5/((c*x**2)**(sympy.S(3)/2)*(a + b*x)**2)
    F = -a**2*x/(b**3*c*sqrt(c*x**2)*(a + b*x)) - 2*a*x*log(a + b*x)/(b**3*c*sqrt(c*x**2)) + x**2/(b**2*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_917():
    f = x**4/((c*x**2)**(sympy.S(3)/2)*(a + b*x)**2)
    F = a*x/(b**2*c*sqrt(c*x**2)*(a + b*x)) + x*log(a + b*x)/(b**2*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_918():
    f = x**3/((c*x**2)**(sympy.S(3)/2)*(a + b*x)**2)
    F = -x/(b*c*sqrt(c*x**2)*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_919():
    f = x**2/((c*x**2)**(sympy.S(3)/2)*(a + b*x)**2)
    F = x/(a*c*sqrt(c*x**2)*(a + b*x)) + x*log(x)/(a**2*c*sqrt(c*x**2)) - x*log(a + b*x)/(a**2*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_920():
    f = x/((c*x**2)**(sympy.S(3)/2)*(a + b*x)**2)
    F = -b*x/(a**2*c*sqrt(c*x**2)*(a + b*x)) - 1/(a**2*c*sqrt(c*x**2)) - 2*b*x*log(x)/(a**3*c*sqrt(c*x**2)) + 2*b*x*log(a + b*x)/(a**3*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_921():
    f = 1/((c*x**2)**(sympy.S(3)/2)*(a + b*x)**2)
    F = -1/(2*a**2*c*x*sqrt(c*x**2)) + b**2*x/(a**3*c*sqrt(c*x**2)*(a + b*x)) + 2*b/(a**3*c*sqrt(c*x**2)) + 3*b**2*x*log(x)/(a**4*c*sqrt(c*x**2)) - 3*b**2*x*log(a + b*x)/(a**4*c*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_922():
    f = x**2*sqrt(c*x**2)*(a + b*x)**n
    F = -a**3*sqrt(c*x**2)*(a + b*x)**(n + 1)/(b**4*x*(n + 1)) + 3*a**2*sqrt(c*x**2)*(a + b*x)**(n + 2)/(b**4*x*(n + 2)) - 3*a*sqrt(c*x**2)*(a + b*x)**(n + 3)/(b**4*x*(n + 3)) + sqrt(c*x**2)*(a + b*x)**(n + 4)/(b**4*x*(n + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_923():
    f = x*sqrt(c*x**2)*(a + b*x)**n
    F = a**2*sqrt(c*x**2)*(a + b*x)**(n + 1)/(b**3*x*(n + 1)) - 2*a*sqrt(c*x**2)*(a + b*x)**(n + 2)/(b**3*x*(n + 2)) + sqrt(c*x**2)*(a + b*x)**(n + 3)/(b**3*x*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_924():
    f = sqrt(c*x**2)*(a + b*x)**n
    F = -a*sqrt(c*x**2)*(a + b*x)**(n + 1)/(b**2*x*(n + 1)) + sqrt(c*x**2)*(a + b*x)**(n + 2)/(b**2*x*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_925():
    f = sqrt(c*x**2)*(a + b*x)**n/x
    F = sqrt(c*x**2)*(a + b*x)**(n + 1)/(b*x*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_926():
    f = sqrt(c*x**2)*(a + b*x)**n/x**2
    F = -sqrt(c*x**2)*(a + b*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*x/a)/(a*x*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_927():
    f = sqrt(c*x**2)*(a + b*x)**n/x**3
    F = b*sqrt(c*x**2)*(a + b*x)**(n + 1)*hyper((2, n + 1), (n + 2,), 1 + b*x/a)/(a**2*x*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_928():
    f = sqrt(c*x**2)*(a + b*x)**n/x**4
    F = -b**2*sqrt(c*x**2)*(a + b*x)**(n + 1)*hyper((3, n + 1), (n + 2,), 1 + b*x/a)/(a**3*x*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_929():
    f = x*(c*x**2)**(sympy.S(3)/2)*(a + b*x)**n
    F = a**4*c*sqrt(c*x**2)*(a + b*x)**(n + 1)/(b**5*x*(n + 1)) - 4*a**3*c*sqrt(c*x**2)*(a + b*x)**(n + 2)/(b**5*x*(n + 2)) + 6*a**2*c*sqrt(c*x**2)*(a + b*x)**(n + 3)/(b**5*x*(n + 3)) - 4*a*c*sqrt(c*x**2)*(a + b*x)**(n + 4)/(b**5*x*(n + 4)) + c*sqrt(c*x**2)*(a + b*x)**(n + 5)/(b**5*x*(n + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_930():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)**n
    F = -a**3*c*sqrt(c*x**2)*(a + b*x)**(n + 1)/(b**4*x*(n + 1)) + 3*a**2*c*sqrt(c*x**2)*(a + b*x)**(n + 2)/(b**4*x*(n + 2)) - 3*a*c*sqrt(c*x**2)*(a + b*x)**(n + 3)/(b**4*x*(n + 3)) + c*sqrt(c*x**2)*(a + b*x)**(n + 4)/(b**4*x*(n + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_931():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)**n/x
    F = a**2*c*sqrt(c*x**2)*(a + b*x)**(n + 1)/(b**3*x*(n + 1)) - 2*a*c*sqrt(c*x**2)*(a + b*x)**(n + 2)/(b**3*x*(n + 2)) + c*sqrt(c*x**2)*(a + b*x)**(n + 3)/(b**3*x*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_932():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)**n/x**2
    F = -a*c*sqrt(c*x**2)*(a + b*x)**(n + 1)/(b**2*x*(n + 1)) + c*sqrt(c*x**2)*(a + b*x)**(n + 2)/(b**2*x*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_933():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)**n/x**3
    F = c*sqrt(c*x**2)*(a + b*x)**(n + 1)/(b*x*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_934():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)**n/x**4
    F = -c*sqrt(c*x**2)*(a + b*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*x/a)/(a*x*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_935():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)**n/x**5
    F = b*c*sqrt(c*x**2)*(a + b*x)**(n + 1)*hyper((2, n + 1), (n + 2,), 1 + b*x/a)/(a**2*x*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_936():
    f = (c*x**2)**(sympy.S(3)/2)*(a + b*x)**n/x**6
    F = -b**2*c*sqrt(c*x**2)*(a + b*x)**(n + 1)*hyper((3, n + 1), (n + 2,), 1 + b*x/a)/(a**3*x*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_937():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)**n
    F = -a**5*c**2*sqrt(c*x**2)*(a + b*x)**(n + 1)/(b**6*x*(n + 1)) + 5*a**4*c**2*sqrt(c*x**2)*(a + b*x)**(n + 2)/(b**6*x*(n + 2)) - 10*a**3*c**2*sqrt(c*x**2)*(a + b*x)**(n + 3)/(b**6*x*(n + 3)) + 10*a**2*c**2*sqrt(c*x**2)*(a + b*x)**(n + 4)/(b**6*x*(n + 4)) - 5*a*c**2*sqrt(c*x**2)*(a + b*x)**(n + 5)/(b**6*x*(n + 5)) + c**2*sqrt(c*x**2)*(a + b*x)**(n + 6)/(b**6*x*(n + 6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_938():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)**n/x
    F = a**4*c**2*sqrt(c*x**2)*(a + b*x)**(n + 1)/(b**5*x*(n + 1)) - 4*a**3*c**2*sqrt(c*x**2)*(a + b*x)**(n + 2)/(b**5*x*(n + 2)) + 6*a**2*c**2*sqrt(c*x**2)*(a + b*x)**(n + 3)/(b**5*x*(n + 3)) - 4*a*c**2*sqrt(c*x**2)*(a + b*x)**(n + 4)/(b**5*x*(n + 4)) + c**2*sqrt(c*x**2)*(a + b*x)**(n + 5)/(b**5*x*(n + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_939():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)**n/x**2
    F = -a**3*c**2*sqrt(c*x**2)*(a + b*x)**(n + 1)/(b**4*x*(n + 1)) + 3*a**2*c**2*sqrt(c*x**2)*(a + b*x)**(n + 2)/(b**4*x*(n + 2)) - 3*a*c**2*sqrt(c*x**2)*(a + b*x)**(n + 3)/(b**4*x*(n + 3)) + c**2*sqrt(c*x**2)*(a + b*x)**(n + 4)/(b**4*x*(n + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_940():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)**n/x**3
    F = a**2*c**2*sqrt(c*x**2)*(a + b*x)**(n + 1)/(b**3*x*(n + 1)) - 2*a*c**2*sqrt(c*x**2)*(a + b*x)**(n + 2)/(b**3*x*(n + 2)) + c**2*sqrt(c*x**2)*(a + b*x)**(n + 3)/(b**3*x*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_941():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)**n/x**4
    F = -a*c**2*sqrt(c*x**2)*(a + b*x)**(n + 1)/(b**2*x*(n + 1)) + c**2*sqrt(c*x**2)*(a + b*x)**(n + 2)/(b**2*x*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_942():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)**n/x**5
    F = c**2*sqrt(c*x**2)*(a + b*x)**(n + 1)/(b*x*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_943():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)**n/x**6
    F = -c**2*sqrt(c*x**2)*(a + b*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*x/a)/(a*x*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_944():
    f = (c*x**2)**(sympy.S(5)/2)*(a + b*x)**n/x**7
    F = b*c**2*sqrt(c*x**2)*(a + b*x)**(n + 1)*hyper((2, n + 1), (n + 2,), 1 + b*x/a)/(a**2*x*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_945():
    f = x**4*(a + b*x)**n/sqrt(c*x**2)
    F = -a**3*x*(a + b*x)**(n + 1)/(b**4*sqrt(c*x**2)*(n + 1)) + 3*a**2*x*(a + b*x)**(n + 2)/(b**4*sqrt(c*x**2)*(n + 2)) - 3*a*x*(a + b*x)**(n + 3)/(b**4*sqrt(c*x**2)*(n + 3)) + x*(a + b*x)**(n + 4)/(b**4*sqrt(c*x**2)*(n + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_946():
    f = x**3*(a + b*x)**n/sqrt(c*x**2)
    F = a**2*x*(a + b*x)**(n + 1)/(b**3*sqrt(c*x**2)*(n + 1)) - 2*a*x*(a + b*x)**(n + 2)/(b**3*sqrt(c*x**2)*(n + 2)) + x*(a + b*x)**(n + 3)/(b**3*sqrt(c*x**2)*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_947():
    f = x**2*(a + b*x)**n/sqrt(c*x**2)
    F = -a*x*(a + b*x)**(n + 1)/(b**2*sqrt(c*x**2)*(n + 1)) + x*(a + b*x)**(n + 2)/(b**2*sqrt(c*x**2)*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_948():
    f = x*(a + b*x)**n/sqrt(c*x**2)
    F = x*(a + b*x)**(n + 1)/(b*sqrt(c*x**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_949():
    f = (a + b*x)**n/sqrt(c*x**2)
    F = -x*(a + b*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*x/a)/(a*sqrt(c*x**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_950():
    f = (a + b*x)**n/(x*sqrt(c*x**2))
    F = b*x*(a + b*x)**(n + 1)*hyper((2, n + 1), (n + 2,), 1 + b*x/a)/(a**2*sqrt(c*x**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_951():
    f = (a + b*x)**n/(x**2*sqrt(c*x**2))
    F = -b**2*x*(a + b*x)**(n + 1)*hyper((3, n + 1), (n + 2,), 1 + b*x/a)/(a**3*sqrt(c*x**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_952():
    f = x**6*(a + b*x)**n/(c*x**2)**(sympy.S(3)/2)
    F = -a**3*x*(a + b*x)**(n + 1)/(b**4*c*sqrt(c*x**2)*(n + 1)) + 3*a**2*x*(a + b*x)**(n + 2)/(b**4*c*sqrt(c*x**2)*(n + 2)) - 3*a*x*(a + b*x)**(n + 3)/(b**4*c*sqrt(c*x**2)*(n + 3)) + x*(a + b*x)**(n + 4)/(b**4*c*sqrt(c*x**2)*(n + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_953():
    f = x**5*(a + b*x)**n/(c*x**2)**(sympy.S(3)/2)
    F = a**2*x*(a + b*x)**(n + 1)/(b**3*c*sqrt(c*x**2)*(n + 1)) - 2*a*x*(a + b*x)**(n + 2)/(b**3*c*sqrt(c*x**2)*(n + 2)) + x*(a + b*x)**(n + 3)/(b**3*c*sqrt(c*x**2)*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_954():
    f = x**4*(a + b*x)**n/(c*x**2)**(sympy.S(3)/2)
    F = -a*x*(a + b*x)**(n + 1)/(b**2*c*sqrt(c*x**2)*(n + 1)) + x*(a + b*x)**(n + 2)/(b**2*c*sqrt(c*x**2)*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_955():
    f = x**3*(a + b*x)**n/(c*x**2)**(sympy.S(3)/2)
    F = x*(a + b*x)**(n + 1)/(b*c*sqrt(c*x**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_956():
    f = x**2*(a + b*x)**n/(c*x**2)**(sympy.S(3)/2)
    F = -x*(a + b*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*x/a)/(a*c*sqrt(c*x**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_957():
    f = x*(a + b*x)**n/(c*x**2)**(sympy.S(3)/2)
    F = b*x*(a + b*x)**(n + 1)*hyper((2, n + 1), (n + 2,), 1 + b*x/a)/(a**2*c*sqrt(c*x**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_958():
    f = (a + b*x)**n/(c*x**2)**(sympy.S(3)/2)
    F = -b**2*x*(a + b*x)**(n + 1)*hyper((3, n + 1), (n + 2,), 1 + b*x/a)/(a**3*c*sqrt(c*x**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_959():
    f = (a + b*x)**n/(x*(c*x**2)**(sympy.S(3)/2))
    F = b**3*x*(a + b*x)**(n + 1)*hyper((4, n + 1), (n + 2,), 1 + b*x/a)/(a**4*c*sqrt(c*x**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_960():
    f = x**8*(a + b*x)**n/(c*x**2)**(sympy.S(5)/2)
    F = -a**3*x*(a + b*x)**(n + 1)/(b**4*c**2*sqrt(c*x**2)*(n + 1)) + 3*a**2*x*(a + b*x)**(n + 2)/(b**4*c**2*sqrt(c*x**2)*(n + 2)) - 3*a*x*(a + b*x)**(n + 3)/(b**4*c**2*sqrt(c*x**2)*(n + 3)) + x*(a + b*x)**(n + 4)/(b**4*c**2*sqrt(c*x**2)*(n + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_961():
    f = x**7*(a + b*x)**n/(c*x**2)**(sympy.S(5)/2)
    F = a**2*x*(a + b*x)**(n + 1)/(b**3*c**2*sqrt(c*x**2)*(n + 1)) - 2*a*x*(a + b*x)**(n + 2)/(b**3*c**2*sqrt(c*x**2)*(n + 2)) + x*(a + b*x)**(n + 3)/(b**3*c**2*sqrt(c*x**2)*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_962():
    f = x**6*(a + b*x)**n/(c*x**2)**(sympy.S(5)/2)
    F = -a*x*(a + b*x)**(n + 1)/(b**2*c**2*sqrt(c*x**2)*(n + 1)) + x*(a + b*x)**(n + 2)/(b**2*c**2*sqrt(c*x**2)*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_963():
    f = x**5*(a + b*x)**n/(c*x**2)**(sympy.S(5)/2)
    F = x*(a + b*x)**(n + 1)/(b*c**2*sqrt(c*x**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_964():
    f = x**4*(a + b*x)**n/(c*x**2)**(sympy.S(5)/2)
    F = -x*(a + b*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*x/a)/(a*c**2*sqrt(c*x**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_965():
    f = x**3*(a + b*x)**n/(c*x**2)**(sympy.S(5)/2)
    F = b*x*(a + b*x)**(n + 1)*hyper((2, n + 1), (n + 2,), 1 + b*x/a)/(a**2*c**2*sqrt(c*x**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_966():
    f = x**2*(a + b*x)**n/(c*x**2)**(sympy.S(5)/2)
    F = -b**2*x*(a + b*x)**(n + 1)*hyper((3, n + 1), (n + 2,), 1 + b*x/a)/(a**3*c**2*sqrt(c*x**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_967():
    f = x*(a + b*x)**n/(c*x**2)**(sympy.S(5)/2)
    F = b**3*x*(a + b*x)**(n + 1)*hyper((4, n + 1), (n + 2,), 1 + b*x/a)/(a**4*c**2*sqrt(c*x**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_968():
    f = (c*x**2)**(sympy.S(5)/2)*(d*x)**m*(a + b*x)
    F = a*c**2*sqrt(c*x**2)*(d*x)**(m + 6)/(d**6*x*(m + 6)) + b*c**2*sqrt(c*x**2)*(d*x)**(m + 7)/(d**7*x*(m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_969():
    f = (c*x**2)**(sympy.S(3)/2)*(d*x)**m*(a + b*x)
    F = a*c*sqrt(c*x**2)*(d*x)**(m + 4)/(d**4*x*(m + 4)) + b*c*sqrt(c*x**2)*(d*x)**(m + 5)/(d**5*x*(m + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_970():
    f = sqrt(c*x**2)*(d*x)**m*(a + b*x)
    F = a*sqrt(c*x**2)*(d*x)**(m + 2)/(d**2*x*(m + 2)) + b*sqrt(c*x**2)*(d*x)**(m + 3)/(d**3*x*(m + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_971():
    f = (d*x)**m*(a + b*x)/sqrt(c*x**2)
    F = a*x*(d*x)**m/(m*sqrt(c*x**2)) + b*x*(d*x)**(m + 1)/(d*sqrt(c*x**2)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_972():
    f = (d*x)**m*(a + b*x)/(c*x**2)**(sympy.S(3)/2)
    F = -a*d**2*x*(d*x)**(m - 2)/(c*sqrt(c*x**2)*(2 - m)) - b*d*x*(d*x)**(m - 1)/(c*sqrt(c*x**2)*(1 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_973():
    f = (d*x)**m*(a + b*x)/(c*x**2)**(sympy.S(5)/2)
    F = -a*d**4*x*(d*x)**(m - 4)/(c**2*sqrt(c*x**2)*(4 - m)) - b*d**3*x*(d*x)**(m - 3)/(c**2*sqrt(c*x**2)*(3 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_974():
    f = (c*x**2)**(sympy.S(5)/2)*(d*x)**m*(a + b*x)**2
    F = a**2*c**2*sqrt(c*x**2)*(d*x)**(m + 6)/(d**6*x*(m + 6)) + 2*a*b*c**2*sqrt(c*x**2)*(d*x)**(m + 7)/(d**7*x*(m + 7)) + b**2*c**2*sqrt(c*x**2)*(d*x)**(m + 8)/(d**8*x*(m + 8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_975():
    f = (c*x**2)**(sympy.S(3)/2)*(d*x)**m*(a + b*x)**2
    F = a**2*c*sqrt(c*x**2)*(d*x)**(m + 4)/(d**4*x*(m + 4)) + 2*a*b*c*sqrt(c*x**2)*(d*x)**(m + 5)/(d**5*x*(m + 5)) + b**2*c*sqrt(c*x**2)*(d*x)**(m + 6)/(d**6*x*(m + 6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_976():
    f = sqrt(c*x**2)*(d*x)**m*(a + b*x)**2
    F = a**2*sqrt(c*x**2)*(d*x)**(m + 2)/(d**2*x*(m + 2)) + 2*a*b*sqrt(c*x**2)*(d*x)**(m + 3)/(d**3*x*(m + 3)) + b**2*sqrt(c*x**2)*(d*x)**(m + 4)/(d**4*x*(m + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_977():
    f = (d*x)**m*(a + b*x)**2/sqrt(c*x**2)
    F = a**2*x*(d*x)**m/(m*sqrt(c*x**2)) + 2*a*b*x*(d*x)**(m + 1)/(d*sqrt(c*x**2)*(m + 1)) + b**2*x*(d*x)**(m + 2)/(d**2*sqrt(c*x**2)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_978():
    f = (d*x)**m*(a + b*x)**2/(c*x**2)**(sympy.S(3)/2)
    F = -a**2*d**2*x*(d*x)**(m - 2)/(c*sqrt(c*x**2)*(2 - m)) - 2*a*b*d*x*(d*x)**(m - 1)/(c*sqrt(c*x**2)*(1 - m)) + b**2*x*(d*x)**m/(c*m*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_979():
    f = (d*x)**m*(a + b*x)**2/(c*x**2)**(sympy.S(5)/2)
    F = -a**2*d**4*x*(d*x)**(m - 4)/(c**2*sqrt(c*x**2)*(4 - m)) - 2*a*b*d**3*x*(d*x)**(m - 3)/(c**2*sqrt(c*x**2)*(3 - m)) - b**2*d**2*x*(d*x)**(m - 2)/(c**2*sqrt(c*x**2)*(2 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_980():
    f = (c*x**2)**(sympy.S(5)/2)*(d*x)**m*(a + b*x)**n
    F = c**2*sqrt(c*x**2)*(d*x)**(m + 6)*(a + b*x)**n*hyper((-n, m + 6), (m + 7,), -b*x/a)/(d**6*x*(1 + b*x/a)**n*(m + 6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_981():
    f = (c*x**2)**(sympy.S(3)/2)*(d*x)**m*(a + b*x)**n
    F = c*sqrt(c*x**2)*(d*x)**(m + 4)*(a + b*x)**n*hyper((-n, m + 4), (m + 5,), -b*x/a)/(d**4*x*(1 + b*x/a)**n*(m + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_982():
    f = sqrt(c*x**2)*(d*x)**m*(a + b*x)**n
    F = sqrt(c*x**2)*(d*x)**(m + 2)*(a + b*x)**n*hyper((-n, m + 2), (m + 3,), -b*x/a)/(d**2*x*(1 + b*x/a)**n*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_983():
    f = (d*x)**m*(a + b*x)**n/sqrt(c*x**2)
    F = x*(d*x)**m*(a + b*x)**n*hyper((m, -n), (m + 1,), -b*x/a)/(m*sqrt(c*x**2)*(1 + b*x/a)**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_984():
    f = (d*x)**m*(a + b*x)**n/(c*x**2)**(sympy.S(3)/2)
    F = -d**2*x*(d*x)**(m - 2)*(a + b*x)**n*hyper((-n, m - 2), (m - 1,), -b*x/a)/(c*sqrt(c*x**2)*(1 + b*x/a)**n*(2 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_985():
    f = (d*x)**m*(a + b*x)**n/(c*x**2)**(sympy.S(5)/2)
    F = -d**4*x*(d*x)**(m - 4)*(a + b*x)**n*hyper((-n, m - 4), (m - 3,), -b*x/a)/(c**2*sqrt(c*x**2)*(1 + b*x/a)**n*(4 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_986():
    f = x**3*(c*x**2)**p*(a + b*x)**(-2*p - 5)
    F = x**4*(c*x**2)**p*(a + b*x)**(-2*p - 4)/(2*a*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_987():
    f = x**2*(c*x**2)**p*(a + b*x)**(-2*p - 4)
    F = x**3*(c*x**2)**p*(a + b*x)**(-2*p - 3)/(a*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_988():
    f = x*(c*x**2)**p*(a + b*x)**(-2*p - 3)
    F = x**2*(c*x**2)**p*(a + b*x)**(-2*p - 2)/(2*a*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_989():
    f = (c*x**2)**p*(a + b*x)**(-2*p - 2)
    F = x*(c*x**2)**p*(a + b*x)**(-2*p - 1)/(a*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_990():
    f = (c*x**2)**p*(a + b*x)**(-2*p - 1)/x
    F = (c*x**2)**p/(2*a*p*(a + b*x)**(2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_991():
    f = (c*x**2)**p/(x**2*(a + b*x)**(2*p))
    F = -(c*x**2)**p*(a + b*x)**(1 - 2*p)/(a*x*(1 - 2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_992():
    f = (c*x**2)**p*(a + b*x)**(1 - 2*p)/x**3
    F = -(c*x**2)**p*(a + b*x)**(2 - 2*p)/(2*a*x**2*(1 - p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_993():
    f = (c*x**2)**p*(a + b*x)**(2 - 2*p)/x**4
    F = -(c*x**2)**p*(a + b*x)**(3 - 2*p)/(a*x**3*(3 - 2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_994():
    f = x**m*(c*x**2)**p*(a + b*x)**(-m - 2*p - 2)
    F = x**(m + 1)*(c*x**2)**p*(a + b*x)**(-m - 2*p - 1)/(a*(m + 2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_995():
    f = (c*x**2)**p*(d*x)**m*(a + b*x)**(-m - 2*p - 2)
    F = x*(c*x**2)**p*(d*x)**m*(a + b*x)**(-m - 2*p - 1)/(a*(m + 2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_996():
    f = x**m*(c*x**2)**p*(a + b*x)**n
    F = x**(m + 1)*(c*x**2)**p*(a + b*x)**n*hyper((-n, m + 2*p + 1), (m + 2*p + 2,), -b*x/a)/((1 + b*x/a)**n*(m + 2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_997():
    f = (a + b*x)**5/(a*d/b + d*x)**3
    F = b**2*(a + b*x)**3/(3*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_998():
    f = (a + b*x)**4/(a*d/b + d*x)**3
    F = a*b**3*x/d**3 + b**4*x**2/(2*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_999():
    f = (a + b*x)**3/(a*d/b + d*x)**3
    F = b**3*x/d**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1000():
    f = (a + b*x)**2/(a*d/b + d*x)**3
    F = b**2*log(a + b*x)/d**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1001():
    f = (a + b*x)/(a*d/b + d*x)**3
    F = -b**2/(d**3*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1002():
    f = 1/((a + b*x)*(a*d/b + d*x)**3)
    F = -b**2/(3*d**3*(a + b*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1003():
    f = 1/((a + b*x)**2*(a*d/b + d*x)**3)
    F = -b**2/(4*d**3*(a + b*x)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1004():
    f = 1/((a + b*x)**3*(a*d/b + d*x)**3)
    F = -b**2/(5*d**3*(a + b*x)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1005():
    f = (b*c/d + b*x)**5/(c + d*x)**3
    F = b**5*(c + d*x)**3/(3*d**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1006():
    f = (b*c/d + b*x)**4/(c + d*x)**3
    F = b**4*c*x/d**4 + b**4*x**2/(2*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1007():
    f = (b*c/d + b*x)**3/(c + d*x)**3
    F = b**3*x/d**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1008():
    f = (b*c/d + b*x)**2/(c + d*x)**3
    F = b**2*log(c + d*x)/d**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1009():
    f = (b*c/d + b*x)/(c + d*x)**3
    F = -b/(d**2*(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1010():
    f = 1/((c + d*x)**3*(b*c/d + b*x))
    F = -1/(3*b*(c + d*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1011():
    f = 1/((c + d*x)**3*(b*c/d + b*x)**2)
    F = -d/(4*b**2*(c + d*x)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1012():
    f = 1/((c + d*x)**3*(b*c/d + b*x)**3)
    F = -d**2/(5*b**3*(c + d*x)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1013():
    f = (a + b*x)**5*(a*c + b*c*x)**n
    F = (a*c + b*c*x)**(n + 6)/(b*c**6*(n + 6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1014():
    f = (a + b*x)**5*(a*c + b*c*x)**3
    F = c**3*(a + b*x)**9/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1015():
    f = (a + b*x)**5*(a*c + b*c*x)**2
    F = c**2*(a + b*x)**8/(8*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1016():
    f = (a + b*x)**5*(a*c + b*c*x)
    F = c*(a + b*x)**7/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1017():
    f = (a + b*x)**5/(a*c + b*c*x)
    F = (a + b*x)**5/(5*b*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1018():
    f = (a + b*x)**5/(a*c + b*c*x)**2
    F = (a + b*x)**4/(4*b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1019():
    f = (a + b*x)**5/(a*c + b*c*x)**3
    F = (a + b*x)**3/(3*b*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1020():
    f = (a + b*x)**5/(a*c + b*c*x)**4
    F = a*x/c**4 + b*x**2/(2*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1021():
    f = (a + b*x)**5/(a*c + b*c*x)**5
    F = x/c**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1022():
    f = (a + b*x)**5/(a*c + b*c*x)**6
    F = log(a + b*x)/(b*c**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1023():
    f = (a + b*x)**5/(a*c + b*c*x)**7
    F = -1/(b*c**7*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1024():
    f = (a + b*x)**5/(a*c + b*c*x)**8
    F = -1/(2*b*c**8*(a + b*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1025():
    f = 1/(sqrt(-3*x - 2)*sqrt(3*x + 2))
    F = sqrt(3*x + 2)*log(3*x + 2)/(3*sqrt(-3*x - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1026():
    f = (a + b*x)*(a*c - b*c*x)**3
    F = -a*c**3*(a - b*x)**4/(2*b) + c**3*(a - b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1027():
    f = (a + b*x)*(a*c - b*c*x)**2
    F = -2*a*c**2*(a - b*x)**3/(3*b) + c**2*(a - b*x)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1028():
    f = (a + b*x)*(a*c - b*c*x)
    F = a**2*c*x - b**2*c*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1029():
    f = a + b*x
    F = a*x + b*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1030():
    f = (a + b*x)/(a*c - b*c*x)
    F = -2*a*log(a - b*x)/(b*c) - x/c
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1031():
    f = (a + b*x)/(a*c - b*c*x)**2
    F = 2*a/(b*c**2*(a - b*x)) + log(a - b*x)/(b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1032():
    f = (a + b*x)/(a*c - b*c*x)**3
    F = x/(c**3*(a - b*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1033():
    f = (a + b*x)/(a*c - b*c*x)**4
    F = 2*a/(3*b*c**4*(a - b*x)**3) - 1/(2*b*c**4*(a - b*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1034():
    f = (a + b*x)/(a*c - b*c*x)**5
    F = a/(2*b*c**5*(a - b*x)**4) - 1/(3*b*c**5*(a - b*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1035():
    f = (a + b*x)/(a*c - b*c*x)**6
    F = 2*a/(5*b*c**6*(a - b*x)**5) - 1/(4*b*c**6*(a - b*x)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1036():
    f = (a + b*x)**2*(a*c - b*c*x)**3
    F = -a**2*c**3*(a - b*x)**4/b + 4*a*c**3*(a - b*x)**5/(5*b) - c**3*(a - b*x)**6/(6*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1037():
    f = (a + b*x)**2*(a*c - b*c*x)**2
    F = a**4*c**2*x - 2*a**2*b**2*c**2*x**3/3 + b**4*c**2*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1038():
    f = (a + b*x)**2*(a*c - b*c*x)
    F = 2*a*c*(a + b*x)**3/(3*b) - c*(a + b*x)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1039():
    f = (a + b*x)**2
    F = (a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1040():
    f = (a + b*x)**2/(a*c - b*c*x)
    F = -4*a**2*log(a - b*x)/(b*c) - 2*a*x/c - (a + b*x)**2/(2*b*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1041():
    f = (a + b*x)**2/(a*c - b*c*x)**2
    F = 4*a**2/(b*c**2*(a - b*x)) + 4*a*log(a - b*x)/(b*c**2) + x/c**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1042():
    f = (a + b*x)**2/(a*c - b*c*x)**3
    F = 2*a**2/(b*c**3*(a - b*x)**2) - 4*a/(b*c**3*(a - b*x)) - log(a - b*x)/(b*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1043():
    f = (a + b*x)**2/(a*c - b*c*x)**4
    F = (a + b*x)**3/(6*a*b*c**4*(a - b*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1044():
    f = (a + b*x)**2/(a*c - b*c*x)**5
    F = a**2/(b*c**5*(a - b*x)**4) - 4*a/(3*b*c**5*(a - b*x)**3) + 1/(2*b*c**5*(a - b*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1045():
    f = (a + b*x)**2/(a*c - b*c*x)**6
    F = 4*a**2/(5*b*c**6*(a - b*x)**5) - a/(b*c**6*(a - b*x)**4) + 1/(3*b*c**6*(a - b*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1046():
    f = (a + b*x)**2/(a*c - b*c*x)**7
    F = 2*a**2/(3*b*c**7*(a - b*x)**6) - 4*a/(5*b*c**7*(a - b*x)**5) + 1/(4*b*c**7*(a - b*x)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1047():
    f = (a*c - b*c*x)**3/(a + b*x)
    F = 8*a**3*c**3*log(a + b*x)/b - 4*a**2*c**3*x + a*c**3*(a - b*x)**2/b + c**3*(a - b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1048():
    f = (a*c - b*c*x)**2/(a + b*x)
    F = 4*a**2*c**2*log(a + b*x)/b - 2*a*c**2*x + c**2*(a - b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1049():
    f = (a*c - b*c*x)/(a + b*x)
    F = 2*a*c*log(a + b*x)/b - c*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1050():
    f = 1/(a + b*x)
    F = log(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1051():
    f = 1/((a + b*x)*(a*c - b*c*x))
    F = atanh(b*x/a)/(a*b*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1052():
    f = 1/((a + b*x)*(a*c - b*c*x)**2)
    F = 1/(2*a*b*c**2*(a - b*x)) + atanh(b*x/a)/(2*a**2*b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1053():
    f = 1/((a + b*x)*(a*c - b*c*x)**3)
    F = 1/(4*a*b*c**3*(a - b*x)**2) + 1/(4*a**2*b*c**3*(a - b*x)) + atanh(b*x/a)/(4*a**3*b*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1054():
    f = (a*c - b*c*x)**3/(a + b*x)**2
    F = -8*a**3*c**3/(b*(a + b*x)) - 12*a**2*c**3*log(a + b*x)/b + 5*a*c**3*x - b*c**3*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1055():
    f = (a*c - b*c*x)**2/(a + b*x)**2
    F = -4*a**2*c**2/(b*(a + b*x)) - 4*a*c**2*log(a + b*x)/b + c**2*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1056():
    f = (a*c - b*c*x)/(a + b*x)**2
    F = -2*a*c/(b*(a + b*x)) - c*log(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1057():
    f = (a + b*x)**(-2)
    F = -1/(b*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1058():
    f = 1/((a + b*x)**2*(a*c - b*c*x))
    F = -1/(2*a*b*c*(a + b*x)) + atanh(b*x/a)/(2*a**2*b*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1059():
    f = 1/((a + b*x)**2*(a*c - b*c*x)**2)
    F = x/(2*a**2*c**2*(a**2 - b**2*x**2)) + atanh(b*x/a)/(2*a**3*b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1060():
    f = 1/((a + b*x)**2*(a*c - b*c*x)**3)
    F = 1/(8*a**2*b*c**3*(a - b*x)**2) - 1/(8*a**3*b*c**3*(a + b*x)) + 1/(4*a**3*b*c**3*(a - b*x)) + 3*atanh(b*x/a)/(8*a**4*b*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1061():
    f = (1 - x)**(sympy.S(9)/2)*sqrt(x + 1)
    F = 21*x*sqrt(1 - x)*sqrt(x + 1)/16 + (1 - x)**(sympy.S(9)/2)*(x + 1)**(sympy.S(3)/2)/6 + 3*(1 - x)**(sympy.S(7)/2)*(x + 1)**(sympy.S(3)/2)/10 + 21*(1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(3)/2)/40 + 7*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/8 + 21*asin(x)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1062():
    f = (1 - x)**(sympy.S(7)/2)*sqrt(x + 1)
    F = 7*x*sqrt(1 - x)*sqrt(x + 1)/8 + (1 - x)**(sympy.S(7)/2)*(x + 1)**(sympy.S(3)/2)/5 + 7*(1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(3)/2)/20 + 7*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/12 + 7*asin(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1063():
    f = (1 - x)**(sympy.S(5)/2)*sqrt(x + 1)
    F = 5*x*sqrt(1 - x)*sqrt(x + 1)/8 + (1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(3)/2)/4 + 5*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/12 + 5*asin(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1064():
    f = (1 - x)**(sympy.S(3)/2)*sqrt(x + 1)
    F = x*sqrt(1 - x)*sqrt(x + 1)/2 + (1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/3 + asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1065():
    f = sqrt(1 - x)*sqrt(x + 1)
    F = x*sqrt(1 - x)*sqrt(x + 1)/2 + asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1066():
    f = sqrt(x + 1)/sqrt(1 - x)
    F = -sqrt(1 - x)*sqrt(x + 1) + asin(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1067():
    f = sqrt(x + 1)/(1 - x)**(sympy.S(3)/2)
    F = -asin(x) + 2*sqrt(x + 1)/sqrt(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1068():
    f = sqrt(x + 1)/(1 - x)**(sympy.S(5)/2)
    F = (x + 1)**(sympy.S(3)/2)/(3*(1 - x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1069():
    f = sqrt(x + 1)/(1 - x)**(sympy.S(7)/2)
    F = (x + 1)**(sympy.S(3)/2)/(15*(1 - x)**(sympy.S(3)/2)) + (x + 1)**(sympy.S(3)/2)/(5*(1 - x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1070():
    f = sqrt(x + 1)/(1 - x)**(sympy.S(9)/2)
    F = 2*(x + 1)**(sympy.S(3)/2)/(105*(1 - x)**(sympy.S(3)/2)) + 2*(x + 1)**(sympy.S(3)/2)/(35*(1 - x)**(sympy.S(5)/2)) + (x + 1)**(sympy.S(3)/2)/(7*(1 - x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1071():
    f = sqrt(x + 1)/(1 - x)**(sympy.S(11)/2)
    F = 2*(x + 1)**(sympy.S(3)/2)/(315*(1 - x)**(sympy.S(3)/2)) + 2*(x + 1)**(sympy.S(3)/2)/(105*(1 - x)**(sympy.S(5)/2)) + (x + 1)**(sympy.S(3)/2)/(21*(1 - x)**(sympy.S(7)/2)) + (x + 1)**(sympy.S(3)/2)/(9*(1 - x)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1072():
    f = sqrt(x + 1)/(1 - x)**(sympy.S(13)/2)
    F = 8*(x + 1)**(sympy.S(3)/2)/(3465*(1 - x)**(sympy.S(3)/2)) + 8*(x + 1)**(sympy.S(3)/2)/(1155*(1 - x)**(sympy.S(5)/2)) + 4*(x + 1)**(sympy.S(3)/2)/(231*(1 - x)**(sympy.S(7)/2)) + 4*(x + 1)**(sympy.S(3)/2)/(99*(1 - x)**(sympy.S(9)/2)) + (x + 1)**(sympy.S(3)/2)/(11*(1 - x)**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1073():
    f = (1 - x)**(sympy.S(9)/2)*(x + 1)**(sympy.S(3)/2)
    F = 3*x*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/8 + 9*x*sqrt(1 - x)*sqrt(x + 1)/16 + (1 - x)**(sympy.S(9)/2)*(x + 1)**(sympy.S(5)/2)/7 + 3*(1 - x)**(sympy.S(7)/2)*(x + 1)**(sympy.S(5)/2)/14 + 3*(1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(5)/2)/10 + 9*asin(x)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1074():
    f = (1 - x)**(sympy.S(7)/2)*(x + 1)**(sympy.S(3)/2)
    F = 7*x*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/24 + 7*x*sqrt(1 - x)*sqrt(x + 1)/16 + (1 - x)**(sympy.S(7)/2)*(x + 1)**(sympy.S(5)/2)/6 + 7*(1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(5)/2)/30 + 7*asin(x)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1075():
    f = (1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(3)/2)
    F = x*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/4 + 3*x*sqrt(1 - x)*sqrt(x + 1)/8 + (1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(5)/2)/5 + 3*asin(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1076():
    f = (1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)
    F = x*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/4 + 3*x*sqrt(1 - x)*sqrt(x + 1)/8 + 3*asin(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1077():
    f = sqrt(1 - x)*(x + 1)**(sympy.S(3)/2)
    F = x*sqrt(1 - x)*sqrt(x + 1)/2 - (1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/3 + asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1078():
    f = (x + 1)**(sympy.S(3)/2)/sqrt(1 - x)
    F = -sqrt(1 - x)*(x + 1)**(sympy.S(3)/2)/2 - 3*sqrt(1 - x)*sqrt(x + 1)/2 + 3*asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1079():
    f = (x + 1)**(sympy.S(3)/2)/(1 - x)**(sympy.S(3)/2)
    F = 3*sqrt(1 - x)*sqrt(x + 1) - 3*asin(x) + 2*(x + 1)**(sympy.S(3)/2)/sqrt(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1080():
    f = (x + 1)**(sympy.S(3)/2)/(1 - x)**(sympy.S(5)/2)
    F = asin(x) - 2*sqrt(x + 1)/sqrt(1 - x) + 2*(x + 1)**(sympy.S(3)/2)/(3*(1 - x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1081():
    f = (x + 1)**(sympy.S(3)/2)/(1 - x)**(sympy.S(7)/2)
    F = (x + 1)**(sympy.S(5)/2)/(5*(1 - x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1082():
    f = (x + 1)**(sympy.S(3)/2)/(1 - x)**(sympy.S(9)/2)
    F = (x + 1)**(sympy.S(5)/2)/(35*(1 - x)**(sympy.S(5)/2)) + (x + 1)**(sympy.S(5)/2)/(7*(1 - x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1083():
    f = (x + 1)**(sympy.S(3)/2)/(1 - x)**(sympy.S(11)/2)
    F = 2*(x + 1)**(sympy.S(5)/2)/(315*(1 - x)**(sympy.S(5)/2)) + 2*(x + 1)**(sympy.S(5)/2)/(63*(1 - x)**(sympy.S(7)/2)) + (x + 1)**(sympy.S(5)/2)/(9*(1 - x)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1084():
    f = (x + 1)**(sympy.S(3)/2)/(1 - x)**(sympy.S(13)/2)
    F = 2*(x + 1)**(sympy.S(5)/2)/(1155*(1 - x)**(sympy.S(5)/2)) + 2*(x + 1)**(sympy.S(5)/2)/(231*(1 - x)**(sympy.S(7)/2)) + (x + 1)**(sympy.S(5)/2)/(33*(1 - x)**(sympy.S(9)/2)) + (x + 1)**(sympy.S(5)/2)/(11*(1 - x)**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1085():
    f = (x + 1)**(sympy.S(3)/2)/(1 - x)**(sympy.S(15)/2)
    F = 8*(x + 1)**(sympy.S(5)/2)/(15015*(1 - x)**(sympy.S(5)/2)) + 8*(x + 1)**(sympy.S(5)/2)/(3003*(1 - x)**(sympy.S(7)/2)) + 4*(x + 1)**(sympy.S(5)/2)/(429*(1 - x)**(sympy.S(9)/2)) + 4*(x + 1)**(sympy.S(5)/2)/(143*(1 - x)**(sympy.S(11)/2)) + (x + 1)**(sympy.S(5)/2)/(13*(1 - x)**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1086():
    f = (1 - x)**(sympy.S(11)/2)*(x + 1)**(sympy.S(5)/2)
    F = 11*x*(1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(5)/2)/48 + 55*x*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/192 + 55*x*sqrt(1 - x)*sqrt(x + 1)/128 + (1 - x)**(sympy.S(11)/2)*(x + 1)**(sympy.S(7)/2)/9 + 11*(1 - x)**(sympy.S(9)/2)*(x + 1)**(sympy.S(7)/2)/72 + 11*(1 - x)**(sympy.S(7)/2)*(x + 1)**(sympy.S(7)/2)/56 + 55*asin(x)/128
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1087():
    f = (1 - x)**(sympy.S(9)/2)*(x + 1)**(sympy.S(5)/2)
    F = 3*x*(1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(5)/2)/16 + 15*x*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/64 + 45*x*sqrt(1 - x)*sqrt(x + 1)/128 + (1 - x)**(sympy.S(9)/2)*(x + 1)**(sympy.S(7)/2)/8 + 9*(1 - x)**(sympy.S(7)/2)*(x + 1)**(sympy.S(7)/2)/56 + 45*asin(x)/128
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1088():
    f = (1 - x)**(sympy.S(7)/2)*(x + 1)**(sympy.S(5)/2)
    F = x*(1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(5)/2)/6 + 5*x*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/24 + 5*x*sqrt(1 - x)*sqrt(x + 1)/16 + (1 - x)**(sympy.S(7)/2)*(x + 1)**(sympy.S(7)/2)/7 + 5*asin(x)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1089():
    f = (1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(5)/2)
    F = x*(1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(5)/2)/6 + 5*x*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/24 + 5*x*sqrt(1 - x)*sqrt(x + 1)/16 + 5*asin(x)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1090():
    f = (1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(5)/2)
    F = x*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/4 + 3*x*sqrt(1 - x)*sqrt(x + 1)/8 - (1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(5)/2)/5 + 3*asin(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1091():
    f = sqrt(1 - x)*(x + 1)**(sympy.S(5)/2)
    F = 5*x*sqrt(1 - x)*sqrt(x + 1)/8 - (1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(5)/2)/4 - 5*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/12 + 5*asin(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1092():
    f = (x + 1)**(sympy.S(5)/2)/sqrt(1 - x)
    F = -sqrt(1 - x)*(x + 1)**(sympy.S(5)/2)/3 - 5*sqrt(1 - x)*(x + 1)**(sympy.S(3)/2)/6 - 5*sqrt(1 - x)*sqrt(x + 1)/2 + 5*asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1093():
    f = (x + 1)**(sympy.S(5)/2)/(1 - x)**(sympy.S(3)/2)
    F = 5*sqrt(1 - x)*(x + 1)**(sympy.S(3)/2)/2 + 15*sqrt(1 - x)*sqrt(x + 1)/2 - 15*asin(x)/2 + 2*(x + 1)**(sympy.S(5)/2)/sqrt(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1094():
    f = (x + 1)**(sympy.S(5)/2)/(1 - x)**(sympy.S(5)/2)
    F = -5*sqrt(1 - x)*sqrt(x + 1) + 5*asin(x) - 10*(x + 1)**(sympy.S(3)/2)/(3*sqrt(1 - x)) + 2*(x + 1)**(sympy.S(5)/2)/(3*(1 - x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1095():
    f = (x + 1)**(sympy.S(5)/2)/(1 - x)**(sympy.S(7)/2)
    F = -asin(x) + 2*sqrt(x + 1)/sqrt(1 - x) - 2*(x + 1)**(sympy.S(3)/2)/(3*(1 - x)**(sympy.S(3)/2)) + 2*(x + 1)**(sympy.S(5)/2)/(5*(1 - x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1096():
    f = (x + 1)**(sympy.S(5)/2)/(1 - x)**(sympy.S(9)/2)
    F = (x + 1)**(sympy.S(7)/2)/(7*(1 - x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1097():
    f = (x + 1)**(sympy.S(5)/2)/(1 - x)**(sympy.S(11)/2)
    F = (x + 1)**(sympy.S(7)/2)/(63*(1 - x)**(sympy.S(7)/2)) + (x + 1)**(sympy.S(7)/2)/(9*(1 - x)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1098():
    f = (x + 1)**(sympy.S(5)/2)/(1 - x)**(sympy.S(13)/2)
    F = 2*(x + 1)**(sympy.S(7)/2)/(693*(1 - x)**(sympy.S(7)/2)) + 2*(x + 1)**(sympy.S(7)/2)/(99*(1 - x)**(sympy.S(9)/2)) + (x + 1)**(sympy.S(7)/2)/(11*(1 - x)**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1099():
    f = (x + 1)**(sympy.S(5)/2)/(1 - x)**(sympy.S(15)/2)
    F = 2*(x + 1)**(sympy.S(7)/2)/(3003*(1 - x)**(sympy.S(7)/2)) + 2*(x + 1)**(sympy.S(7)/2)/(429*(1 - x)**(sympy.S(9)/2)) + 3*(x + 1)**(sympy.S(7)/2)/(143*(1 - x)**(sympy.S(11)/2)) + (x + 1)**(sympy.S(7)/2)/(13*(1 - x)**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1100():
    f = (x + 1)**(sympy.S(5)/2)/(1 - x)**(sympy.S(17)/2)
    F = 8*(x + 1)**(sympy.S(7)/2)/(45045*(1 - x)**(sympy.S(7)/2)) + 8*(x + 1)**(sympy.S(7)/2)/(6435*(1 - x)**(sympy.S(9)/2)) + 4*(x + 1)**(sympy.S(7)/2)/(715*(1 - x)**(sympy.S(11)/2)) + 4*(x + 1)**(sympy.S(7)/2)/(195*(1 - x)**(sympy.S(13)/2)) + (x + 1)**(sympy.S(7)/2)/(15*(1 - x)**(sympy.S(15)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1101():
    f = (x + 1)**(sympy.S(5)/2)/(1 - x)**(sympy.S(19)/2)
    F = 8*(x + 1)**(sympy.S(7)/2)/(153153*(1 - x)**(sympy.S(7)/2)) + 8*(x + 1)**(sympy.S(7)/2)/(21879*(1 - x)**(sympy.S(9)/2)) + 4*(x + 1)**(sympy.S(7)/2)/(2431*(1 - x)**(sympy.S(11)/2)) + 4*(x + 1)**(sympy.S(7)/2)/(663*(1 - x)**(sympy.S(13)/2)) + (x + 1)**(sympy.S(7)/2)/(51*(1 - x)**(sympy.S(15)/2)) + (x + 1)**(sympy.S(7)/2)/(17*(1 - x)**(sympy.S(17)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1102():
    f = (a*x + 1)**(sympy.S(3)/2)/sqrt(-a*x + 1)
    F = -sqrt(-a*x + 1)*(a*x + 1)**(sympy.S(3)/2)/(2*a) - 3*sqrt(-a*x + 1)*sqrt(a*x + 1)/(2*a) + 3*asin(a*x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1103():
    f = (a*x + 1)*sqrt(-a**2*x**2 + 1)/(-a*x + 1)
    F = -3*sqrt(-a**2*x**2 + 1)/(2*a) + 3*asin(a*x)/(2*a) - (-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1104():
    f = (1 - x)**(sympy.S(7)/2)/sqrt(x + 1)
    F = (1 - x)**(sympy.S(7)/2)*sqrt(x + 1)/4 + 7*(1 - x)**(sympy.S(5)/2)*sqrt(x + 1)/12 + 35*(1 - x)**(sympy.S(3)/2)*sqrt(x + 1)/24 + 35*sqrt(1 - x)*sqrt(x + 1)/8 + 35*asin(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1105():
    f = (1 - x)**(sympy.S(5)/2)/sqrt(x + 1)
    F = (1 - x)**(sympy.S(5)/2)*sqrt(x + 1)/3 + 5*(1 - x)**(sympy.S(3)/2)*sqrt(x + 1)/6 + 5*sqrt(1 - x)*sqrt(x + 1)/2 + 5*asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1106():
    f = (1 - x)**(sympy.S(3)/2)/sqrt(x + 1)
    F = (1 - x)**(sympy.S(3)/2)*sqrt(x + 1)/2 + 3*sqrt(1 - x)*sqrt(x + 1)/2 + 3*asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1107():
    f = sqrt(1 - x)/sqrt(x + 1)
    F = sqrt(1 - x)*sqrt(x + 1) + asin(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1108():
    f = 1/(sqrt(1 - x)*sqrt(x + 1))
    F = asin(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1109():
    f = 1/((1 - x)**(sympy.S(3)/2)*sqrt(x + 1))
    F = sqrt(x + 1)/sqrt(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1110():
    f = 1/((1 - x)**(sympy.S(5)/2)*sqrt(x + 1))
    F = sqrt(x + 1)/(3*sqrt(1 - x)) + sqrt(x + 1)/(3*(1 - x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1111():
    f = 1/((1 - x)**(sympy.S(7)/2)*sqrt(x + 1))
    F = 2*sqrt(x + 1)/(15*sqrt(1 - x)) + 2*sqrt(x + 1)/(15*(1 - x)**(sympy.S(3)/2)) + sqrt(x + 1)/(5*(1 - x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1112():
    f = 1/((1 - x)**(sympy.S(9)/2)*sqrt(x + 1))
    F = 2*sqrt(x + 1)/(35*sqrt(1 - x)) + 2*sqrt(x + 1)/(35*(1 - x)**(sympy.S(3)/2)) + 3*sqrt(x + 1)/(35*(1 - x)**(sympy.S(5)/2)) + sqrt(x + 1)/(7*(1 - x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1113():
    f = 1/((1 - x)**(sympy.S(11)/2)*sqrt(x + 1))
    F = 8*sqrt(x + 1)/(315*sqrt(1 - x)) + 8*sqrt(x + 1)/(315*(1 - x)**(sympy.S(3)/2)) + 4*sqrt(x + 1)/(105*(1 - x)**(sympy.S(5)/2)) + 4*sqrt(x + 1)/(63*(1 - x)**(sympy.S(7)/2)) + sqrt(x + 1)/(9*(1 - x)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1114():
    f = (1 - x)**(sympy.S(7)/2)/(x + 1)**(sympy.S(3)/2)
    F = -2*(1 - x)**(sympy.S(7)/2)/sqrt(x + 1) - 7*(1 - x)**(sympy.S(5)/2)*sqrt(x + 1)/3 - 35*(1 - x)**(sympy.S(3)/2)*sqrt(x + 1)/6 - 35*sqrt(1 - x)*sqrt(x + 1)/2 - 35*asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1115():
    f = (1 - x)**(sympy.S(5)/2)/(x + 1)**(sympy.S(3)/2)
    F = -2*(1 - x)**(sympy.S(5)/2)/sqrt(x + 1) - 5*(1 - x)**(sympy.S(3)/2)*sqrt(x + 1)/2 - 15*sqrt(1 - x)*sqrt(x + 1)/2 - 15*asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1116():
    f = (1 - x)**(sympy.S(3)/2)/(x + 1)**(sympy.S(3)/2)
    F = -2*(1 - x)**(sympy.S(3)/2)/sqrt(x + 1) - 3*sqrt(1 - x)*sqrt(x + 1) - 3*asin(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1117():
    f = sqrt(1 - x)/(x + 1)**(sympy.S(3)/2)
    F = -2*sqrt(1 - x)/sqrt(x + 1) - asin(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1118():
    f = 1/(sqrt(1 - x)*(x + 1)**(sympy.S(3)/2))
    F = -sqrt(1 - x)/sqrt(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1119():
    f = 1/((1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2))
    F = x/(sqrt(1 - x)*sqrt(x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1120():
    f = 1/((1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(3)/2))
    F = 2*x/(3*sqrt(1 - x)*sqrt(x + 1)) + 1/(3*(1 - x)**(sympy.S(3)/2)*sqrt(x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1121():
    f = 1/((1 - x)**(sympy.S(7)/2)*(x + 1)**(sympy.S(3)/2))
    F = 2*x/(5*sqrt(1 - x)*sqrt(x + 1)) + 1/(5*(1 - x)**(sympy.S(3)/2)*sqrt(x + 1)) + 1/(5*(1 - x)**(sympy.S(5)/2)*sqrt(x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1122():
    f = 1/((1 - x)**(sympy.S(9)/2)*(x + 1)**(sympy.S(3)/2))
    F = 8*x/(35*sqrt(1 - x)*sqrt(x + 1)) + 4/(35*(1 - x)**(sympy.S(3)/2)*sqrt(x + 1)) + 4/(35*(1 - x)**(sympy.S(5)/2)*sqrt(x + 1)) + 1/(7*(1 - x)**(sympy.S(7)/2)*sqrt(x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1123():
    f = 1/((1 - x)**(sympy.S(11)/2)*(x + 1)**(sympy.S(3)/2))
    F = 8*x/(63*sqrt(1 - x)*sqrt(x + 1)) + 4/(63*(1 - x)**(sympy.S(3)/2)*sqrt(x + 1)) + 4/(63*(1 - x)**(sympy.S(5)/2)*sqrt(x + 1)) + 5/(63*(1 - x)**(sympy.S(7)/2)*sqrt(x + 1)) + 1/(9*(1 - x)**(sympy.S(9)/2)*sqrt(x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1124():
    f = (1 - x)**(sympy.S(9)/2)/(x + 1)**(sympy.S(5)/2)
    F = -2*(1 - x)**(sympy.S(9)/2)/(3*(x + 1)**(sympy.S(3)/2)) + 6*(1 - x)**(sympy.S(7)/2)/sqrt(x + 1) + 7*(1 - x)**(sympy.S(5)/2)*sqrt(x + 1) + 35*(1 - x)**(sympy.S(3)/2)*sqrt(x + 1)/2 + 105*sqrt(1 - x)*sqrt(x + 1)/2 + 105*asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1125():
    f = (1 - x)**(sympy.S(7)/2)/(x + 1)**(sympy.S(5)/2)
    F = -2*(1 - x)**(sympy.S(7)/2)/(3*(x + 1)**(sympy.S(3)/2)) + 14*(1 - x)**(sympy.S(5)/2)/(3*sqrt(x + 1)) + 35*(1 - x)**(sympy.S(3)/2)*sqrt(x + 1)/6 + 35*sqrt(1 - x)*sqrt(x + 1)/2 + 35*asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1126():
    f = (1 - x)**(sympy.S(5)/2)/(x + 1)**(sympy.S(5)/2)
    F = -2*(1 - x)**(sympy.S(5)/2)/(3*(x + 1)**(sympy.S(3)/2)) + 10*(1 - x)**(sympy.S(3)/2)/(3*sqrt(x + 1)) + 5*sqrt(1 - x)*sqrt(x + 1) + 5*asin(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1127():
    f = (1 - x)**(sympy.S(3)/2)/(x + 1)**(sympy.S(5)/2)
    F = -2*(1 - x)**(sympy.S(3)/2)/(3*(x + 1)**(sympy.S(3)/2)) + 2*sqrt(1 - x)/sqrt(x + 1) + asin(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1128():
    f = sqrt(1 - x)/(x + 1)**(sympy.S(5)/2)
    F = -(1 - x)**(sympy.S(3)/2)/(3*(x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1129():
    f = 1/(sqrt(1 - x)*(x + 1)**(sympy.S(5)/2))
    F = -sqrt(1 - x)/(3*sqrt(x + 1)) - sqrt(1 - x)/(3*(x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1130():
    f = 1/((1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(5)/2))
    F = -2*sqrt(1 - x)/(3*sqrt(x + 1)) - 2*sqrt(1 - x)/(3*(x + 1)**(sympy.S(3)/2)) + 1/(sqrt(1 - x)*(x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1131():
    f = 1/((1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(5)/2))
    F = 2*x/(3*sqrt(1 - x)*sqrt(x + 1)) + x/(3*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1132():
    f = 1/((1 - x)**(sympy.S(7)/2)*(x + 1)**(sympy.S(5)/2))
    F = 8*x/(15*sqrt(1 - x)*sqrt(x + 1)) + 4*x/(15*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)) + 1/(5*(1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1133():
    f = 1/((1 - x)**(sympy.S(9)/2)*(x + 1)**(sympy.S(5)/2))
    F = 8*x/(21*sqrt(1 - x)*sqrt(x + 1)) + 4*x/(21*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)) + 1/(7*(1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(3)/2)) + 1/(7*(1 - x)**(sympy.S(7)/2)*(x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1134():
    f = 1/((1 - x)**(sympy.S(11)/2)*(x + 1)**(sympy.S(5)/2))
    F = 16*x/(63*sqrt(1 - x)*sqrt(x + 1)) + 8*x/(63*(1 - x)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)) + 2/(21*(1 - x)**(sympy.S(5)/2)*(x + 1)**(sympy.S(3)/2)) + 2/(21*(1 - x)**(sympy.S(7)/2)*(x + 1)**(sympy.S(3)/2)) + 1/(9*(1 - x)**(sympy.S(9)/2)*(x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1135():
    f = (a*x + a)**(sympy.S(5)/2)*(-c*x + c)**(sympy.S(5)/2)
    F = 5*a**(sympy.S(5)/2)*c**(sympy.S(5)/2)*atan(sqrt(c)*sqrt(a*x + a)/(sqrt(a)*sqrt(-c*x + c)))/8 + 5*a**2*c**2*x*sqrt(a*x + a)*sqrt(-c*x + c)/16 + 5*a*c*x*(a*x + a)**(sympy.S(3)/2)*(-c*x + c)**(sympy.S(3)/2)/24 + x*(a*x + a)**(sympy.S(5)/2)*(-c*x + c)**(sympy.S(5)/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1136():
    f = (a*x + a)**(sympy.S(3)/2)*(-c*x + c)**(sympy.S(3)/2)
    F = 3*a**(sympy.S(3)/2)*c**(sympy.S(3)/2)*atan(sqrt(c)*sqrt(a*x + a)/(sqrt(a)*sqrt(-c*x + c)))/4 + 3*a*c*x*sqrt(a*x + a)*sqrt(-c*x + c)/8 + x*(a*x + a)**(sympy.S(3)/2)*(-c*x + c)**(sympy.S(3)/2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1137():
    f = sqrt(a*x + a)*sqrt(-c*x + c)
    F = sqrt(a)*sqrt(c)*atan(sqrt(c)*sqrt(a*x + a)/(sqrt(a)*sqrt(-c*x + c))) + x*sqrt(a*x + a)*sqrt(-c*x + c)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1138():
    f = 1/(sqrt(a*x + a)*sqrt(-c*x + c))
    F = 2*atan(sqrt(c)*sqrt(a*x + a)/(sqrt(a)*sqrt(-c*x + c)))/(sqrt(a)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1139():
    f = 1/((a*x + a)**(sympy.S(3)/2)*(-c*x + c)**(sympy.S(3)/2))
    F = x/(a*c*sqrt(a*x + a)*sqrt(-c*x + c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1140():
    f = 1/((a*x + a)**(sympy.S(5)/2)*(-c*x + c)**(sympy.S(5)/2))
    F = x/(3*a*c*(a*x + a)**(sympy.S(3)/2)*(-c*x + c)**(sympy.S(3)/2)) + 2*x/(3*a**2*c**2*sqrt(a*x + a)*sqrt(-c*x + c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1141():
    f = 1/((a*x + a)**(sympy.S(7)/2)*(-c*x + c)**(sympy.S(7)/2))
    F = x/(5*a*c*(a*x + a)**(sympy.S(5)/2)*(-c*x + c)**(sympy.S(5)/2)) + 4*x/(15*a**2*c**2*(a*x + a)**(sympy.S(3)/2)*(-c*x + c)**(sympy.S(3)/2)) + 8*x/(15*a**3*c**3*sqrt(a*x + a)*sqrt(-c*x + c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1142():
    f = 1/((a*x + a)**(sympy.S(9)/2)*(-c*x + c)**(sympy.S(9)/2))
    F = x/(7*a*c*(a*x + a)**(sympy.S(7)/2)*(-c*x + c)**(sympy.S(7)/2)) + 6*x/(35*a**2*c**2*(a*x + a)**(sympy.S(5)/2)*(-c*x + c)**(sympy.S(5)/2)) + 8*x/(35*a**3*c**3*(a*x + a)**(sympy.S(3)/2)*(-c*x + c)**(sympy.S(3)/2)) + 16*x/(35*a**4*c**4*sqrt(a*x + a)*sqrt(-c*x + c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1143():
    f = (a + b*x)**(sympy.S(5)/2)*(a*c - b*c*x)**(sympy.S(5)/2)
    F = 5*a**6*c**(sympy.S(5)/2)*atan(sqrt(c)*sqrt(a + b*x)/sqrt(c*(a - b*x)))/(8*b) + 5*a**4*c**2*x*sqrt(a + b*x)*sqrt(a*c - b*c*x)/16 + 5*a**2*c*x*(a + b*x)**(sympy.S(3)/2)*(a*c - b*c*x)**(sympy.S(3)/2)/24 + x*(a + b*x)**(sympy.S(5)/2)*(a*c - b*c*x)**(sympy.S(5)/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1144():
    f = (a + b*x)**(sympy.S(3)/2)*(a*c - b*c*x)**(sympy.S(3)/2)
    F = 3*a**4*c**(sympy.S(3)/2)*atan(sqrt(c)*sqrt(a + b*x)/sqrt(c*(a - b*x)))/(4*b) + 3*a**2*c*x*sqrt(a + b*x)*sqrt(a*c - b*c*x)/8 + x*(a + b*x)**(sympy.S(3)/2)*(a*c - b*c*x)**(sympy.S(3)/2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1145():
    f = sqrt(a + b*x)*sqrt(a*c - b*c*x)
    F = a**2*sqrt(c)*atan(sqrt(c)*sqrt(a + b*x)/sqrt(c*(a - b*x)))/b + x*sqrt(a + b*x)*sqrt(a*c - b*c*x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1146():
    f = 1/(sqrt(a + b*x)*sqrt(a*c - b*c*x))
    F = 2*atan(sqrt(c)*sqrt(a + b*x)/sqrt(c*(a - b*x)))/(b*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1147():
    f = 1/((a + b*x)**(sympy.S(3)/2)*(a*c - b*c*x)**(sympy.S(3)/2))
    F = x/(a**2*c*sqrt(a + b*x)*sqrt(a*c - b*c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1148():
    f = 1/((a + b*x)**(sympy.S(5)/2)*(a*c - b*c*x)**(sympy.S(5)/2))
    F = x/(3*a**2*c*(a + b*x)**(sympy.S(3)/2)*(a*c - b*c*x)**(sympy.S(3)/2)) + 2*x/(3*a**4*c**2*sqrt(a + b*x)*sqrt(a*c - b*c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1149():
    f = 1/((a + b*x)**(sympy.S(7)/2)*(a*c - b*c*x)**(sympy.S(7)/2))
    F = x/(5*a**2*c*(a + b*x)**(sympy.S(5)/2)*(a*c - b*c*x)**(sympy.S(5)/2)) + 4*x/(15*a**4*c**2*(a + b*x)**(sympy.S(3)/2)*(a*c - b*c*x)**(sympy.S(3)/2)) + 8*x/(15*a**6*c**3*sqrt(a + b*x)*sqrt(a*c - b*c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1150():
    f = 1/((a + b*x)**(sympy.S(9)/2)*(a*c - b*c*x)**(sympy.S(9)/2))
    F = x/(7*a**2*c*(a + b*x)**(sympy.S(7)/2)*(a*c - b*c*x)**(sympy.S(7)/2)) + 6*x/(35*a**4*c**2*(a + b*x)**(sympy.S(5)/2)*(a*c - b*c*x)**(sympy.S(5)/2)) + 8*x/(35*a**6*c**3*(a + b*x)**(sympy.S(3)/2)*(a*c - b*c*x)**(sympy.S(3)/2)) + 16*x/(35*a**8*c**4*sqrt(a + b*x)*sqrt(a*c - b*c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1151():
    f = (3 - 6*x)**(sympy.S(5)/2)*(4*x + 2)**(sympy.S(5)/2)
    F = 6*sqrt(6)*x*(1 - 2*x)**(sympy.S(5)/2)*(2*x + 1)**(sympy.S(5)/2) + 15*sqrt(6)*x*(1 - 2*x)**(sympy.S(3)/2)*(2*x + 1)**(sympy.S(3)/2)/2 + 45*sqrt(6)*x*sqrt(1 - 2*x)*sqrt(2*x + 1)/4 + 45*sqrt(6)*asin(2*x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1152():
    f = (3 - 6*x)**(sympy.S(3)/2)*(4*x + 2)**(sympy.S(3)/2)
    F = 3*sqrt(6)*x*(1 - 2*x)**(sympy.S(3)/2)*(2*x + 1)**(sympy.S(3)/2)/2 + 9*sqrt(6)*x*sqrt(1 - 2*x)*sqrt(2*x + 1)/4 + 9*sqrt(6)*asin(2*x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1153():
    f = sqrt(3 - 6*x)*sqrt(4*x + 2)
    F = sqrt(6)*x*sqrt(1 - 2*x)*sqrt(2*x + 1)/2 + sqrt(6)*asin(2*x)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1154():
    f = 1/(sqrt(3 - 6*x)*sqrt(4*x + 2))
    F = sqrt(6)*asin(2*x)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1155():
    f = 1/((3 - 6*x)**(sympy.S(3)/2)*(4*x + 2)**(sympy.S(3)/2))
    F = sqrt(6)*x/(36*sqrt(1 - 2*x)*sqrt(2*x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1156():
    f = 1/((3 - 6*x)**(sympy.S(5)/2)*(4*x + 2)**(sympy.S(5)/2))
    F = sqrt(6)*x/(324*sqrt(1 - 2*x)*sqrt(2*x + 1)) + sqrt(6)*x/(648*(1 - 2*x)**(sympy.S(3)/2)*(2*x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1157():
    f = 1/((3 - 6*x)**(sympy.S(7)/2)*(4*x + 2)**(sympy.S(7)/2))
    F = sqrt(6)*x/(2430*sqrt(1 - 2*x)*sqrt(2*x + 1)) + sqrt(6)*x/(4860*(1 - 2*x)**(sympy.S(3)/2)*(2*x + 1)**(sympy.S(3)/2)) + sqrt(6)*x/(6480*(1 - 2*x)**(sympy.S(5)/2)*(2*x + 1)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1158():
    f = (3 - x)**(sympy.S(3)/2)*(x - 2)**(sympy.S(3)/2)
    F = -(3 - x)**(sympy.S(5)/2)*(x - 2)**(sympy.S(3)/2)/4 - (3 - x)**(sympy.S(5)/2)*sqrt(x - 2)/8 + (3 - x)**(sympy.S(3)/2)*sqrt(x - 2)/32 + 3*sqrt(3 - x)*sqrt(x - 2)/64 + 3*asin(2*x - 5)/128
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1159():
    f = sqrt(3 - x)*sqrt(x - 2)
    F = -(3 - x)**(sympy.S(3)/2)*sqrt(x - 2)/2 + sqrt(3 - x)*sqrt(x - 2)/4 + asin(2*x - 5)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1160():
    f = 1/(sqrt(3 - x)*sqrt(x - 2))
    F = asin(2*x - 5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1161():
    f = 1/((3 - x)**(sympy.S(3)/2)*(x - 2)**(sympy.S(3)/2))
    F = -4*sqrt(3 - x)/sqrt(x - 2) + 2/(sqrt(3 - x)*sqrt(x - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1162():
    f = 1/((3 - x)**(sympy.S(5)/2)*(x - 2)**(sympy.S(5)/2))
    F = -32*sqrt(3 - x)/(3*sqrt(x - 2)) - 16*sqrt(3 - x)/(3*(x - 2)**(sympy.S(3)/2)) + 4/(sqrt(3 - x)*(x - 2)**(sympy.S(3)/2)) + 2/(3*(3 - x)**(sympy.S(3)/2)*(x - 2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1163():
    f = 1/((3 - x)**(sympy.S(3)/2)*(x + 3)**(sympy.S(3)/2))
    F = x/(9*sqrt(3 - x)*sqrt(x + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1164():
    f = 1/((-b*x + 3)**(sympy.S(3)/2)*(b*x + 3)**(sympy.S(3)/2))
    F = x/(9*sqrt(-b*x + 3)*sqrt(b*x + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1165():
    f = 1/((6 - 2*x)**(sympy.S(3)/2)*(x + 3)**(sympy.S(3)/2))
    F = sqrt(2)*x/(36*sqrt(3 - x)*sqrt(x + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1166():
    f = 1/((-2*b*x + 6)**(sympy.S(3)/2)*(b*x + 3)**(sympy.S(3)/2))
    F = sqrt(2)*x/(36*sqrt(-b*x + 3)*sqrt(b*x + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1167():
    f = 1/(sqrt(a + b*x)*sqrt(-a*d + b*d*x))
    F = 2*atanh(sqrt(d)*sqrt(a + b*x)/sqrt(-a*d + b*d*x))/(b*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1168():
    f = 1/((-3*e*x + 6)**(sympy.S(1)/4)*(e*x + 2)**(sympy.S(3)/4))
    F = -sqrt(2)*3**(sympy.S(3)/4)*log((sqrt(-3*e*x + 6) - sqrt(6)*(-e*x + 2)**(sympy.S(1)/4)*(e*x + 2)**(sympy.S(1)/4) + sqrt(3)*sqrt(e*x + 2))/sqrt(e*x + 2))/(6*e) + sqrt(2)*3**(sympy.S(3)/4)*log((sqrt(-3*e*x + 6) + sqrt(6)*(-e*x + 2)**(sympy.S(1)/4)*(e*x + 2)**(sympy.S(1)/4) + sqrt(3)*sqrt(e*x + 2))/sqrt(e*x + 2))/(6*e) - sqrt(2)*3**(sympy.S(3)/4)*atan(sqrt(2)*(-e*x + 2)**(sympy.S(1)/4)/(e*x + 2)**(sympy.S(1)/4) - 1)/(3*e) - sqrt(2)*3**(sympy.S(3)/4)*atan(sqrt(2)*(-e*x + 2)**(sympy.S(1)/4)/(e*x + 2)**(sympy.S(1)/4) + 1)/(3*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1169():
    f = (-I*a*x + a)**(sympy.S(7)/4)/(I*a*x + a)**(sympy.S(1)/4)
    F = 14*a**2*x/(5*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) - 14*a**2*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(5*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) - 14*I*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4)/15 - 2*I*(-I*a*x + a)**(sympy.S(7)/4)*(I*a*x + a)**(sympy.S(3)/4)/(5*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1170():
    f = (-I*a*x + a)**(sympy.S(3)/4)/(I*a*x + a)**(sympy.S(1)/4)
    F = 2*a*x/((-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) - 2*a*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/((-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) - 2*I*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4)/(3*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1171():
    f = 1/((-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    F = 2*x/((-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) - 2*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/((-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1172():
    f = 1/((-I*a*x + a)**(sympy.S(5)/4)*(I*a*x + a)**(sympy.S(1)/4))
    F = 2*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(a*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) - 2*I/(a*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1173():
    f = 1/((-I*a*x + a)**(sympy.S(9)/4)*(I*a*x + a)**(sympy.S(1)/4))
    F = -4*I/(5*a*(-I*a*x + a)**(sympy.S(5)/4)*(I*a*x + a)**(sympy.S(1)/4)) + 2*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(5*a**2*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1174():
    f = 1/((-I*a*x + a)**(sympy.S(13)/4)*(I*a*x + a)**(sympy.S(1)/4))
    F = -4*I/(15*a**2*(-I*a*x + a)**(sympy.S(5)/4)*(I*a*x + a)**(sympy.S(1)/4)) - 2*I*(I*a*x + a)**(sympy.S(3)/4)/(9*a**2*(-I*a*x + a)**(sympy.S(9)/4)) + 2*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(15*a**3*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1175():
    f = 1/((-I*a*x + a)**(sympy.S(17)/4)*(I*a*x + a)**(sympy.S(1)/4))
    F = -2*I*(I*a*x + a)**(sympy.S(3)/4)/(13*a**2*(-I*a*x + a)**(sympy.S(13)/4)) - 4*I/(39*a**3*(-I*a*x + a)**(sympy.S(5)/4)*(I*a*x + a)**(sympy.S(1)/4)) - 10*I*(I*a*x + a)**(sympy.S(3)/4)/(117*a**3*(-I*a*x + a)**(sympy.S(9)/4)) + 2*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(39*a**4*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1176():
    f = (-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4)
    F = -sqrt(2)*I*log(-sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/4 + sqrt(2)*I*log(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/4 + sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) - 1)/2 + sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + 1)/2 - I*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(3)/4)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1177():
    f = 1/((-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(1)/4))
    F = -sqrt(2)*I*log(-sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/(2*a) + sqrt(2)*I*log(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/(2*a) + sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) - 1)/a + sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + 1)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1178():
    f = 1/((-I*a*x + a)**(sympy.S(7)/4)*(I*a*x + a)**(sympy.S(1)/4))
    F = -2*I*(I*a*x + a)**(sympy.S(3)/4)/(3*a**2*(-I*a*x + a)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1179():
    f = 1/((-I*a*x + a)**(sympy.S(11)/4)*(I*a*x + a)**(sympy.S(1)/4))
    F = -2*I*(I*a*x + a)**(sympy.S(3)/4)/(7*a**2*(-I*a*x + a)**(sympy.S(7)/4)) - 4*I*(I*a*x + a)**(sympy.S(3)/4)/(21*a**3*(-I*a*x + a)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1180():
    f = 1/((-I*a*x + a)**(sympy.S(15)/4)*(I*a*x + a)**(sympy.S(1)/4))
    F = -2*I*(I*a*x + a)**(sympy.S(3)/4)/(11*a**2*(-I*a*x + a)**(sympy.S(11)/4)) - 8*I*(I*a*x + a)**(sympy.S(3)/4)/(77*a**3*(-I*a*x + a)**(sympy.S(7)/4)) - 16*I*(I*a*x + a)**(sympy.S(3)/4)/(231*a**4*(-I*a*x + a)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1181():
    f = 1/((-I*a*x + a)**(sympy.S(19)/4)*(I*a*x + a)**(sympy.S(1)/4))
    F = -2*I*(I*a*x + a)**(sympy.S(3)/4)/(15*a**2*(-I*a*x + a)**(sympy.S(15)/4)) - 4*I*(I*a*x + a)**(sympy.S(3)/4)/(55*a**3*(-I*a*x + a)**(sympy.S(11)/4)) - 16*I*(I*a*x + a)**(sympy.S(3)/4)/(385*a**4*(-I*a*x + a)**(sympy.S(7)/4)) - 32*I*(I*a*x + a)**(sympy.S(3)/4)/(1155*a**5*(-I*a*x + a)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1182():
    f = (-I*a*x + a)**(sympy.S(3)/4)/(I*a*x + a)**(sympy.S(3)/4)
    F = 3*sqrt(2)*I*log(-sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/4 - 3*sqrt(2)*I*log(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/4 + 3*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) - 1)/2 + 3*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + 1)/2 - I*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(1)/4)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1183():
    f = 1/((-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(3)/4))
    F = sqrt(2)*I*log(-sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/(2*a) - sqrt(2)*I*log(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/(2*a) + sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) - 1)/a + sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + 1)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1184():
    f = 1/((-I*a*x + a)**(sympy.S(5)/4)*(I*a*x + a)**(sympy.S(3)/4))
    F = -2*I*(I*a*x + a)**(sympy.S(1)/4)/(a**2*(-I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1185():
    f = 1/((-I*a*x + a)**(sympy.S(9)/4)*(I*a*x + a)**(sympy.S(3)/4))
    F = -2*I*(I*a*x + a)**(sympy.S(1)/4)/(5*a**2*(-I*a*x + a)**(sympy.S(5)/4)) - 4*I*(I*a*x + a)**(sympy.S(1)/4)/(5*a**3*(-I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1186():
    f = 1/((-I*a*x + a)**(sympy.S(13)/4)*(I*a*x + a)**(sympy.S(3)/4))
    F = -2*I*(I*a*x + a)**(sympy.S(1)/4)/(9*a**2*(-I*a*x + a)**(sympy.S(9)/4)) - 8*I*(I*a*x + a)**(sympy.S(1)/4)/(45*a**3*(-I*a*x + a)**(sympy.S(5)/4)) - 16*I*(I*a*x + a)**(sympy.S(1)/4)/(45*a**4*(-I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1187():
    f = (-I*a*x + a)**(sympy.S(5)/4)/(I*a*x + a)**(sympy.S(3)/4)
    F = 10*a**2*(x**2 + 1)**(sympy.S(3)/4)*elliptic_f(atan(x)/2, 2)/(3*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4)) - 10*I*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)/3 - 2*I*(-I*a*x + a)**(sympy.S(5)/4)*(I*a*x + a)**(sympy.S(1)/4)/(3*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1188():
    f = (-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(3)/4)
    F = 2*a*(x**2 + 1)**(sympy.S(3)/4)*elliptic_f(atan(x)/2, 2)/((-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4)) - 2*I*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1189():
    f = 1/((-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4))
    F = 2*(x**2 + 1)**(sympy.S(3)/4)*elliptic_f(atan(x)/2, 2)/((-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1190():
    f = 1/((-I*a*x + a)**(sympy.S(7)/4)*(I*a*x + a)**(sympy.S(3)/4))
    F = 2*(x**2 + 1)**(sympy.S(3)/4)*elliptic_f(atan(x)/2, 2)/(3*a*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4)) - 2*I*(I*a*x + a)**(sympy.S(1)/4)/(3*a**2*(-I*a*x + a)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1191():
    f = 1/((-I*a*x + a)**(sympy.S(11)/4)*(I*a*x + a)**(sympy.S(3)/4))
    F = 2*(x**2 + 1)**(sympy.S(3)/4)*elliptic_f(atan(x)/2, 2)/(7*a**2*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4)) - 2*I*(I*a*x + a)**(sympy.S(1)/4)/(7*a**2*(-I*a*x + a)**(sympy.S(7)/4)) - 2*I*(I*a*x + a)**(sympy.S(1)/4)/(7*a**3*(-I*a*x + a)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1192():
    f = (-I*a*x + a)**(sympy.S(7)/4)/(I*a*x + a)**(sympy.S(7)/4)
    F = -7*sqrt(2)*I*log(-sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/4 + 7*sqrt(2)*I*log(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/4 - 7*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) - 1)/2 - 7*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + 1)/2 + 4*I*(-I*a*x + a)**(sympy.S(7)/4)/(3*a*(I*a*x + a)**(sympy.S(3)/4)) + 7*I*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(1)/4)/(3*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1193():
    f = (-I*a*x + a)**(sympy.S(3)/4)/(I*a*x + a)**(sympy.S(7)/4)
    F = 4*I*(-I*a*x + a)**(sympy.S(3)/4)/(3*a*(I*a*x + a)**(sympy.S(3)/4)) - sqrt(2)*I*log(-sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/(2*a) + sqrt(2)*I*log(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/(2*a) - sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) - 1)/a - sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + 1)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1194():
    f = 1/((-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(7)/4))
    F = 2*I*(-I*a*x + a)**(sympy.S(3)/4)/(3*a**2*(I*a*x + a)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1195():
    f = 1/((-I*a*x + a)**(sympy.S(5)/4)*(I*a*x + a)**(sympy.S(7)/4))
    F = -2*I/(a**2*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(3)/4)) + 4*I*(-I*a*x + a)**(sympy.S(3)/4)/(3*a**3*(I*a*x + a)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1196():
    f = 1/((-I*a*x + a)**(sympy.S(9)/4)*(I*a*x + a)**(sympy.S(7)/4))
    F = -2*I/(5*a**2*(-I*a*x + a)**(sympy.S(5)/4)*(I*a*x + a)**(sympy.S(3)/4)) - 8*I/(5*a**3*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(3)/4)) + 16*I*(-I*a*x + a)**(sympy.S(3)/4)/(15*a**4*(I*a*x + a)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1197():
    f = (-I*a*x + a)**(sympy.S(9)/4)/(I*a*x + a)**(sympy.S(7)/4)
    F = -10*a**2*(x**2 + 1)**(sympy.S(3)/4)*elliptic_f(atan(x)/2, 2)/((-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4)) + 10*I*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4) + 4*I*(-I*a*x + a)**(sympy.S(9)/4)/(3*a*(I*a*x + a)**(sympy.S(3)/4)) + 2*I*(-I*a*x + a)**(sympy.S(5)/4)*(I*a*x + a)**(sympy.S(1)/4)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1198():
    f = (-I*a*x + a)**(sympy.S(5)/4)/(I*a*x + a)**(sympy.S(7)/4)
    F = -10*a*(x**2 + 1)**(sympy.S(3)/4)*elliptic_f(atan(x)/2, 2)/(3*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4)) + 4*I*(-I*a*x + a)**(sympy.S(5)/4)/(3*a*(I*a*x + a)**(sympy.S(3)/4)) + 10*I*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)/(3*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1199():
    f = (-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(7)/4)
    F = -2*(x**2 + 1)**(sympy.S(3)/4)*elliptic_f(atan(x)/2, 2)/(3*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4)) + 4*I*(-I*a*x + a)**(sympy.S(1)/4)/(3*a*(I*a*x + a)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1200():
    f = 1/((-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(7)/4))
    F = 2*(x**2 + 1)**(sympy.S(3)/4)*elliptic_f(atan(x)/2, 2)/(3*a*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4)) + 2*I*(-I*a*x + a)**(sympy.S(1)/4)/(3*a**2*(I*a*x + a)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1201():
    f = 1/((-I*a*x + a)**(sympy.S(7)/4)*(I*a*x + a)**(sympy.S(7)/4))
    F = 2*x/(3*a**2*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4)) + 2*(x**2 + 1)**(sympy.S(3)/4)*elliptic_f(atan(x)/2, 2)/(3*a**2*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1202():
    f = 1/((-I*a*x + a)**(sympy.S(11)/4)*(I*a*x + a)**(sympy.S(7)/4))
    F = -2*I/(7*a**2*(-I*a*x + a)**(sympy.S(7)/4)*(I*a*x + a)**(sympy.S(3)/4)) + 10*x/(21*a**3*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4)) + 10*(x**2 + 1)**(sympy.S(3)/4)*elliptic_f(atan(x)/2, 2)/(21*a**3*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1203():
    f = 1/((-I*a*x + a)**(sympy.S(15)/4)*(I*a*x + a)**(sympy.S(7)/4))
    F = -2*I/(11*a**2*(-I*a*x + a)**(sympy.S(11)/4)*(I*a*x + a)**(sympy.S(3)/4)) - 2*I/(11*a**3*(-I*a*x + a)**(sympy.S(7)/4)*(I*a*x + a)**(sympy.S(3)/4)) + 10*x/(33*a**4*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4)) + 10*(x**2 + 1)**(sympy.S(3)/4)*elliptic_f(atan(x)/2, 2)/(33*a**4*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1204():
    f = (-I*a*x + a)**(sympy.S(7)/4)/(I*a*x + a)**(sympy.S(5)/4)
    F = -14*a*x/((-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) + 14*a*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/((-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) + 4*I*(-I*a*x + a)**(sympy.S(7)/4)/(a*(I*a*x + a)**(sympy.S(1)/4)) + 14*I*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(3)/4)/(3*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1205():
    f = (-I*a*x + a)**(sympy.S(3)/4)/(I*a*x + a)**(sympy.S(5)/4)
    F = -6*x/((-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) + 6*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/((-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) + 4*I*(-I*a*x + a)**(sympy.S(3)/4)/(a*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1206():
    f = 1/((-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(5)/4))
    F = 2*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(a*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) + 2*I/(a*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1207():
    f = 1/((-I*a*x + a)**(sympy.S(5)/4)*(I*a*x + a)**(sympy.S(5)/4))
    F = 2*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(a**2*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1208():
    f = 1/((-I*a*x + a)**(sympy.S(9)/4)*(I*a*x + a)**(sympy.S(5)/4))
    F = -2*I/(5*a**2*(-I*a*x + a)**(sympy.S(5)/4)*(I*a*x + a)**(sympy.S(1)/4)) + 6*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(5*a**3*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1209():
    f = 1/((-I*a*x + a)**(sympy.S(13)/4)*(I*a*x + a)**(sympy.S(5)/4))
    F = -2*I/(9*a**2*(-I*a*x + a)**(sympy.S(9)/4)*(I*a*x + a)**(sympy.S(1)/4)) - 2*I/(9*a**3*(-I*a*x + a)**(sympy.S(5)/4)*(I*a*x + a)**(sympy.S(1)/4)) + 2*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(3*a**4*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1210():
    f = (-I*a*x + a)**(sympy.S(5)/4)/(I*a*x + a)**(sympy.S(5)/4)
    F = 5*sqrt(2)*I*log(-sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/4 - 5*sqrt(2)*I*log(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/4 - 5*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) - 1)/2 - 5*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + 1)/2 + 4*I*(-I*a*x + a)**(sympy.S(5)/4)/(a*(I*a*x + a)**(sympy.S(1)/4)) + 5*I*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(3)/4)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1211():
    f = (-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(5)/4)
    F = 4*I*(-I*a*x + a)**(sympy.S(1)/4)/(a*(I*a*x + a)**(sympy.S(1)/4)) + sqrt(2)*I*log(-sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/(2*a) - sqrt(2)*I*log(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/(2*a) - sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) - 1)/a - sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + 1)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1212():
    f = 1/((-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(5)/4))
    F = 2*I*(-I*a*x + a)**(sympy.S(1)/4)/(a**2*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1213():
    f = 1/((-I*a*x + a)**(sympy.S(7)/4)*(I*a*x + a)**(sympy.S(5)/4))
    F = -2*I/(3*a**2*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(1)/4)) + 4*I*(-I*a*x + a)**(sympy.S(1)/4)/(3*a**3*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1214():
    f = 1/((-I*a*x + a)**(sympy.S(11)/4)*(I*a*x + a)**(sympy.S(5)/4))
    F = -2*I/(7*a**2*(-I*a*x + a)**(sympy.S(7)/4)*(I*a*x + a)**(sympy.S(1)/4)) - 8*I/(21*a**3*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(1)/4)) + 16*I*(-I*a*x + a)**(sympy.S(1)/4)/(21*a**4*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1215():
    f = (-I*a*x + a)**(sympy.S(7)/4)/(I*a*x + a)**(sympy.S(9)/4)
    F = 42*x/(5*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) - 42*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(5*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) + 4*I*(-I*a*x + a)**(sympy.S(7)/4)/(5*a*(I*a*x + a)**(sympy.S(5)/4)) - 28*I*(-I*a*x + a)**(sympy.S(3)/4)/(5*a*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1216():
    f = (-I*a*x + a)**(sympy.S(3)/4)/(I*a*x + a)**(sympy.S(9)/4)
    F = -6*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(5*a*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) + 4*I*(-I*a*x + a)**(sympy.S(3)/4)/(5*a*(I*a*x + a)**(sympy.S(5)/4)) - 6*I/(5*a*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1217():
    f = 1/((-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(9)/4))
    F = 4*I/(5*a*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(5)/4)) + 2*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(5*a**2*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1218():
    f = 1/((-I*a*x + a)**(sympy.S(5)/4)*(I*a*x + a)**(sympy.S(9)/4))
    F = 2*I/(5*a**2*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(5)/4)) + 6*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(5*a**3*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1219():
    f = 1/((-I*a*x + a)**(sympy.S(9)/4)*(I*a*x + a)**(sympy.S(9)/4))
    F = 2*x/(5*a**4*(x**2 + 1)*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) + 6*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(5*a**4*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1220():
    f = 1/((-I*a*x + a)**(sympy.S(13)/4)*(I*a*x + a)**(sympy.S(9)/4))
    F = -2*I/(9*a**2*(-I*a*x + a)**(sympy.S(9)/4)*(I*a*x + a)**(sympy.S(5)/4)) + 14*x/(45*a**5*(x**2 + 1)*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) + 14*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(15*a**5*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1221():
    f = 1/((-I*a*x + a)**(sympy.S(17)/4)*(I*a*x + a)**(sympy.S(9)/4))
    F = -2*I/(13*a**2*(-I*a*x + a)**(sympy.S(13)/4)*(I*a*x + a)**(sympy.S(5)/4)) - 2*I/(13*a**3*(-I*a*x + a)**(sympy.S(9)/4)*(I*a*x + a)**(sympy.S(5)/4)) + 14*x/(65*a**6*(x**2 + 1)*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4)) + 42*(x**2 + 1)**(sympy.S(1)/4)*elliptic_e(atan(x)/2, 2)/(65*a**6*(-I*a*x + a)**(sympy.S(1)/4)*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1222():
    f = (-I*a*x + a)**(sympy.S(5)/4)/(I*a*x + a)**(sympy.S(9)/4)
    F = 4*I*(-I*a*x + a)**(sympy.S(5)/4)/(5*a*(I*a*x + a)**(sympy.S(5)/4)) - 4*I*(-I*a*x + a)**(sympy.S(1)/4)/(a*(I*a*x + a)**(sympy.S(1)/4)) - sqrt(2)*I*log(-sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/(2*a) + sqrt(2)*I*log(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + sqrt(-I*a*x + a)/sqrt(I*a*x + a) + 1)/(2*a) + sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) - 1)/a + sqrt(2)*I*atan(sqrt(2)*(-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(1)/4) + 1)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1223():
    f = (-I*a*x + a)**(sympy.S(1)/4)/(I*a*x + a)**(sympy.S(9)/4)
    F = 2*I*(-I*a*x + a)**(sympy.S(5)/4)/(5*a**2*(I*a*x + a)**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1224():
    f = 1/((-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(9)/4))
    F = 2*I*(-I*a*x + a)**(sympy.S(1)/4)/(5*a**2*(I*a*x + a)**(sympy.S(5)/4)) + 4*I*(-I*a*x + a)**(sympy.S(1)/4)/(5*a**3*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1225():
    f = 1/((-I*a*x + a)**(sympy.S(7)/4)*(I*a*x + a)**(sympy.S(9)/4))
    F = -2*I/(3*a**2*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(5)/4)) + 8*I*(-I*a*x + a)**(sympy.S(1)/4)/(15*a**3*(I*a*x + a)**(sympy.S(5)/4)) + 16*I*(-I*a*x + a)**(sympy.S(1)/4)/(15*a**4*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1226():
    f = 1/((-I*a*x + a)**(sympy.S(11)/4)*(I*a*x + a)**(sympy.S(9)/4))
    F = -2*I/(7*a**2*(-I*a*x + a)**(sympy.S(7)/4)*(I*a*x + a)**(sympy.S(5)/4)) - 4*I/(7*a**3*(-I*a*x + a)**(sympy.S(3)/4)*(I*a*x + a)**(sympy.S(5)/4)) + 16*I*(-I*a*x + a)**(sympy.S(1)/4)/(35*a**4*(I*a*x + a)**(sympy.S(5)/4)) + 32*I*(-I*a*x + a)**(sympy.S(1)/4)/(35*a**5*(I*a*x + a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1227():
    f = (a + b*x)**2*(a*c - b*c*x)**n
    F = -4*a**2*(a*c - b*c*x)**(n + 1)/(b*c*(n + 1)) + 4*a*(a*c - b*c*x)**(n + 2)/(b*c**2*(n + 2)) - (a*c - b*c*x)**(n + 3)/(b*c**3*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1228():
    f = (a + b*x)*(a*c - b*c*x)**n
    F = -2*a*(a*c - b*c*x)**(n + 1)/(b*c*(n + 1)) + (a*c - b*c*x)**(n + 2)/(b*c**2*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1229():
    f = (a*c - b*c*x)**n/(a + b*x)
    F = -(a*c - b*c*x)**(n + 1)*hyper((1, n + 1), (n + 2,), (a - b*x)/(2*a))/(2*a*b*c*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1230():
    f = (a*c - b*c*x)**n/(a + b*x)**2
    F = -(a*c - b*c*x)**(n + 1)*hyper((2, n + 1), (n + 2,), (a - b*x)/(2*a))/(4*a**2*b*c*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1231():
    f = (a*x + a)**m*(-c*x + c)**m
    F = x*(a*x + a)**m*(-c*x + c)**m*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), x**2)/(1 - x**2)**m
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1232():
    f = (a + b*x)**m*(a*c - b*c*x)**m
    F = x*(a + b*x)**m*(a*c - b*c*x)**m*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), b**2*x**2/a**2)/(1 - b**2*x**2/a**2)**m
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1233():
    f = (3 - 6*x)**m*(4*x + 2)**m
    F = 6**m*x*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), 4*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1234():
    f = (a + b*x)**4*(c + d*x)
    F = d*(a + b*x)**6/(6*b**2) + (a + b*x)**5*(-a*d + b*c)/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1235():
    f = (a + b*x)**3*(c + d*x)
    F = d*(a + b*x)**5/(5*b**2) + (a + b*x)**4*(-a*d + b*c)/(4*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1236():
    f = (a + b*x)**2*(c + d*x)
    F = d*(a + b*x)**4/(4*b**2) + (a + b*x)**3*(-a*d + b*c)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1237():
    f = (a + b*x)*(c + d*x)
    F = a*c*x + b*d*x**3/3 + x**2*(a*d/2 + b*c/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1238():
    f = c + d*x
    F = c*x + d*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1239():
    f = (c + d*x)/(a + b*x)
    F = d*x/b + (-a*d + b*c)*log(a + b*x)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1240():
    f = (c + d*x)/(a + b*x)**2
    F = d*log(a + b*x)/b**2 - (-a*d + b*c)/(b**2*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1241():
    f = (c + d*x)/(a + b*x)**3
    F = -(c + d*x)**2/((a + b*x)**2*(-2*a*d + 2*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1242():
    f = (c + d*x)/(a + b*x)**4
    F = -d/(2*b**2*(a + b*x)**2) - (-a*d + b*c)/(3*b**2*(a + b*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1243():
    f = (c + d*x)/(a + b*x)**5
    F = -d/(3*b**2*(a + b*x)**3) - (-a*d + b*c)/(4*b**2*(a + b*x)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1244():
    f = (a + b*x)**4*(c + d*x)**2
    F = d**2*(a + b*x)**7/(7*b**3) + d*(a + b*x)**6*(-a*d + b*c)/(3*b**3) + (a + b*x)**5*(-a*d + b*c)**2/(5*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1245():
    f = (a + b*x)**3*(c + d*x)**2
    F = d**2*(a + b*x)**6/(6*b**3) + 2*d*(a + b*x)**5*(-a*d + b*c)/(5*b**3) + (a + b*x)**4*(-a*d + b*c)**2/(4*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1246():
    f = (a + b*x)**2*(c + d*x)**2
    F = d**2*(a + b*x)**5/(5*b**3) + d*(a + b*x)**4*(-a*d + b*c)/(2*b**3) + (a + b*x)**3*(-a*d + b*c)**2/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1247():
    f = (a + b*x)*(c + d*x)**2
    F = b*(c + d*x)**4/(4*d**2) - (c + d*x)**3*(-a*d + b*c)/(3*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1248():
    f = (c + d*x)**2
    F = (c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1249():
    f = (c + d*x)**2/(a + b*x)
    F = (c + d*x)**2/(2*b) + d*x*(-a*d + b*c)/b**2 + (-a*d + b*c)**2*log(a + b*x)/b**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1250():
    f = (c + d*x)**2/(a + b*x)**2
    F = d**2*x/b**2 + 2*d*(-a*d + b*c)*log(a + b*x)/b**3 - (-a*d + b*c)**2/(b**3*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1251():
    f = (c + d*x)**2/(a + b*x)**3
    F = d**2*log(a + b*x)/b**3 - 2*d*(-a*d + b*c)/(b**3*(a + b*x)) - (-a*d + b*c)**2/(2*b**3*(a + b*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1252():
    f = (c + d*x)**2/(a + b*x)**4
    F = -(c + d*x)**3/((a + b*x)**3*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1253():
    f = (c + d*x)**2/(a + b*x)**5
    F = -d**2/(2*b**3*(a + b*x)**2) - 2*d*(-a*d + b*c)/(3*b**3*(a + b*x)**3) - (-a*d + b*c)**2/(4*b**3*(a + b*x)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1254():
    f = (c + d*x)**2/(a + b*x)**6
    F = -d**2/(3*b**3*(a + b*x)**3) - d*(-a*d + b*c)/(2*b**3*(a + b*x)**4) - (-a*d + b*c)**2/(5*b**3*(a + b*x)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1255():
    f = (c + d*x)**2/(a + b*x)**7
    F = -d**2/(4*b**3*(a + b*x)**4) - 2*d*(-a*d + b*c)/(5*b**3*(a + b*x)**5) - (-a*d + b*c)**2/(6*b**3*(a + b*x)**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1256():
    f = (a + b*x)**5*(c + d*x)**3
    F = d**3*(a + b*x)**9/(9*b**4) + 3*d**2*(a + b*x)**8*(-a*d + b*c)/(8*b**4) + 3*d*(a + b*x)**7*(-a*d + b*c)**2/(7*b**4) + (a + b*x)**6*(-a*d + b*c)**3/(6*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1257():
    f = (a + b*x)**4*(c + d*x)**3
    F = d**3*(a + b*x)**8/(8*b**4) + 3*d**2*(a + b*x)**7*(-a*d + b*c)/(7*b**4) + d*(a + b*x)**6*(-a*d + b*c)**2/(2*b**4) + (a + b*x)**5*(-a*d + b*c)**3/(5*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1258():
    f = (a + b*x)**3*(c + d*x)**3
    F = d**3*(a + b*x)**7/(7*b**4) + d**2*(a + b*x)**6*(-a*d + b*c)/(2*b**4) + 3*d*(a + b*x)**5*(-a*d + b*c)**2/(5*b**4) + (a + b*x)**4*(-a*d + b*c)**3/(4*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1259():
    f = (a + b*x)**2*(c + d*x)**3
    F = b**2*(c + d*x)**6/(6*d**3) - 2*b*(c + d*x)**5*(-a*d + b*c)/(5*d**3) + (c + d*x)**4*(-a*d + b*c)**2/(4*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1260():
    f = (a + b*x)*(c + d*x)**3
    F = b*(c + d*x)**5/(5*d**2) - (c + d*x)**4*(-a*d + b*c)/(4*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1261():
    f = (c + d*x)**3
    F = (c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1262():
    f = (c + d*x)**3/(a + b*x)
    F = (c + d*x)**3/(3*b) + (c + d*x)**2*(-a*d + b*c)/(2*b**2) + d*x*(-a*d + b*c)**2/b**3 + (-a*d + b*c)**3*log(a + b*x)/b**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1263():
    f = (c + d*x)**3/(a + b*x)**2
    F = d**3*x**2/(2*b**2) + d**2*x*(-2*a*d + 3*b*c)/b**3 + 3*d*(-a*d + b*c)**2*log(a + b*x)/b**4 - (-a*d + b*c)**3/(b**4*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1264():
    f = (c + d*x)**3/(a + b*x)**3
    F = d**3*x/b**3 + 3*d**2*(-a*d + b*c)*log(a + b*x)/b**4 - 3*d*(-a*d + b*c)**2/(b**4*(a + b*x)) - (-a*d + b*c)**3/(2*b**4*(a + b*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1265():
    f = (c + d*x)**3/(a + b*x)**4
    F = d**3*log(a + b*x)/b**4 - 3*d**2*(-a*d + b*c)/(b**4*(a + b*x)) - 3*d*(-a*d + b*c)**2/(2*b**4*(a + b*x)**2) - (-a*d + b*c)**3/(3*b**4*(a + b*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1266():
    f = (c + d*x)**3/(a + b*x)**5
    F = -(c + d*x)**4/((a + b*x)**4*(-4*a*d + 4*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1267():
    f = (c + d*x)**3/(a + b*x)**6
    F = d*(c + d*x)**4/(20*(a + b*x)**4*(-a*d + b*c)**2) - (c + d*x)**4/((a + b*x)**5*(-5*a*d + 5*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1268():
    f = (c + d*x)**3/(a + b*x)**7
    F = -d**3/(3*b**4*(a + b*x)**3) - 3*d**2*(-a*d + b*c)/(4*b**4*(a + b*x)**4) - 3*d*(-a*d + b*c)**2/(5*b**4*(a + b*x)**5) - (-a*d + b*c)**3/(6*b**4*(a + b*x)**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1269():
    f = (c + d*x)**3/(a + b*x)**8
    F = -d**3/(4*b**4*(a + b*x)**4) - 3*d**2*(-a*d + b*c)/(5*b**4*(a + b*x)**5) - d*(-a*d + b*c)**2/(2*b**4*(a + b*x)**6) - (-a*d + b*c)**3/(7*b**4*(a + b*x)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1270():
    f = (c + d*x)**3/(a + b*x)**9
    F = -d**3/(5*b**4*(a + b*x)**5) - d**2*(-a*d + b*c)/(2*b**4*(a + b*x)**6) - 3*d*(-a*d + b*c)**2/(7*b**4*(a + b*x)**7) - (-a*d + b*c)**3/(8*b**4*(a + b*x)**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1271():
    f = (a + b*x)**9*(c + d*x)**7
    F = d**7*(a + b*x)**17/(17*b**8) + 7*d**6*(a + b*x)**16*(-a*d + b*c)/(16*b**8) + 7*d**5*(a + b*x)**15*(-a*d + b*c)**2/(5*b**8) + 5*d**4*(a + b*x)**14*(-a*d + b*c)**3/(2*b**8) + 35*d**3*(a + b*x)**13*(-a*d + b*c)**4/(13*b**8) + 7*d**2*(a + b*x)**12*(-a*d + b*c)**5/(4*b**8) + 7*d*(a + b*x)**11*(-a*d + b*c)**6/(11*b**8) + (a + b*x)**10*(-a*d + b*c)**7/(10*b**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1272():
    f = (a + b*x)**8*(c + d*x)**7
    F = d**7*(a + b*x)**16/(16*b**8) + 7*d**6*(a + b*x)**15*(-a*d + b*c)/(15*b**8) + 3*d**5*(a + b*x)**14*(-a*d + b*c)**2/(2*b**8) + 35*d**4*(a + b*x)**13*(-a*d + b*c)**3/(13*b**8) + 35*d**3*(a + b*x)**12*(-a*d + b*c)**4/(12*b**8) + 21*d**2*(a + b*x)**11*(-a*d + b*c)**5/(11*b**8) + 7*d*(a + b*x)**10*(-a*d + b*c)**6/(10*b**8) + (a + b*x)**9*(-a*d + b*c)**7/(9*b**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1273():
    f = (a + b*x)**7*(c + d*x)**7
    F = d**7*(a + b*x)**15/(15*b**8) + d**6*(a + b*x)**14*(-a*d + b*c)/(2*b**8) + 21*d**5*(a + b*x)**13*(-a*d + b*c)**2/(13*b**8) + 35*d**4*(a + b*x)**12*(-a*d + b*c)**3/(12*b**8) + 35*d**3*(a + b*x)**11*(-a*d + b*c)**4/(11*b**8) + 21*d**2*(a + b*x)**10*(-a*d + b*c)**5/(10*b**8) + 7*d*(a + b*x)**9*(-a*d + b*c)**6/(9*b**8) + (a + b*x)**8*(-a*d + b*c)**7/(8*b**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1274():
    f = (a + b*x)**6*(c + d*x)**7
    F = b**6*(c + d*x)**14/(14*d**7) - 6*b**5*(c + d*x)**13*(-a*d + b*c)/(13*d**7) + 5*b**4*(c + d*x)**12*(-a*d + b*c)**2/(4*d**7) - 20*b**3*(c + d*x)**11*(-a*d + b*c)**3/(11*d**7) + 3*b**2*(c + d*x)**10*(-a*d + b*c)**4/(2*d**7) - 2*b*(c + d*x)**9*(-a*d + b*c)**5/(3*d**7) + (c + d*x)**8*(-a*d + b*c)**6/(8*d**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1275():
    f = (a + b*x)**5*(c + d*x)**7
    F = b**5*(c + d*x)**13/(13*d**6) - 5*b**4*(c + d*x)**12*(-a*d + b*c)/(12*d**6) + 10*b**3*(c + d*x)**11*(-a*d + b*c)**2/(11*d**6) - b**2*(c + d*x)**10*(-a*d + b*c)**3/d**6 + 5*b*(c + d*x)**9*(-a*d + b*c)**4/(9*d**6) - (c + d*x)**8*(-a*d + b*c)**5/(8*d**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1276():
    f = (a + b*x)**4*(c + d*x)**7
    F = b**4*(c + d*x)**12/(12*d**5) - 4*b**3*(c + d*x)**11*(-a*d + b*c)/(11*d**5) + 3*b**2*(c + d*x)**10*(-a*d + b*c)**2/(5*d**5) - 4*b*(c + d*x)**9*(-a*d + b*c)**3/(9*d**5) + (c + d*x)**8*(-a*d + b*c)**4/(8*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1277():
    f = (a + b*x)**3*(c + d*x)**7
    F = b**3*(c + d*x)**11/(11*d**4) - 3*b**2*(c + d*x)**10*(-a*d + b*c)/(10*d**4) + b*(c + d*x)**9*(-a*d + b*c)**2/(3*d**4) - (c + d*x)**8*(-a*d + b*c)**3/(8*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1278():
    f = (a + b*x)**2*(c + d*x)**7
    F = b**2*(c + d*x)**10/(10*d**3) - 2*b*(c + d*x)**9*(-a*d + b*c)/(9*d**3) + (c + d*x)**8*(-a*d + b*c)**2/(8*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1279():
    f = (a + b*x)*(c + d*x)**7
    F = b*(c + d*x)**9/(9*d**2) - (c + d*x)**8*(-a*d + b*c)/(8*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1280():
    f = (c + d*x)**7
    F = (c + d*x)**8/(8*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1281():
    f = (c + d*x)**7/(a + b*x)
    F = (c + d*x)**7/(7*b) + (c + d*x)**6*(-a*d + b*c)/(6*b**2) + (c + d*x)**5*(-a*d + b*c)**2/(5*b**3) + (c + d*x)**4*(-a*d + b*c)**3/(4*b**4) + (c + d*x)**3*(-a*d + b*c)**4/(3*b**5) + (c + d*x)**2*(-a*d + b*c)**5/(2*b**6) + d*x*(-a*d + b*c)**6/b**7 + (-a*d + b*c)**7*log(a + b*x)/b**8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1282():
    f = (c + d*x)**7/(a + b*x)**2
    F = 21*d**2*x*(-a*d + b*c)**5/b**7 + d**7*(a + b*x)**6/(6*b**8) + 7*d**6*(a + b*x)**5*(-a*d + b*c)/(5*b**8) + 21*d**5*(a + b*x)**4*(-a*d + b*c)**2/(4*b**8) + 35*d**4*(a + b*x)**3*(-a*d + b*c)**3/(3*b**8) + 35*d**3*(a + b*x)**2*(-a*d + b*c)**4/(2*b**8) + 7*d*(-a*d + b*c)**6*log(a + b*x)/b**8 - (-a*d + b*c)**7/(b**8*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1283():
    f = (c + d*x)**7/(a + b*x)**3
    F = 35*d**3*x*(-a*d + b*c)**4/b**7 + d**7*(a + b*x)**5/(5*b**8) + 7*d**6*(a + b*x)**4*(-a*d + b*c)/(4*b**8) + 7*d**5*(a + b*x)**3*(-a*d + b*c)**2/b**8 + 35*d**4*(a + b*x)**2*(-a*d + b*c)**3/(2*b**8) + 21*d**2*(-a*d + b*c)**5*log(a + b*x)/b**8 - 7*d*(-a*d + b*c)**6/(b**8*(a + b*x)) - (-a*d + b*c)**7/(2*b**8*(a + b*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1284():
    f = (c + d*x)**7/(a + b*x)**4
    F = 35*d**4*x*(-a*d + b*c)**3/b**7 + d**7*(a + b*x)**4/(4*b**8) + 7*d**6*(a + b*x)**3*(-a*d + b*c)/(3*b**8) + 21*d**5*(a + b*x)**2*(-a*d + b*c)**2/(2*b**8) + 35*d**3*(-a*d + b*c)**4*log(a + b*x)/b**8 - 21*d**2*(-a*d + b*c)**5/(b**8*(a + b*x)) - 7*d*(-a*d + b*c)**6/(2*b**8*(a + b*x)**2) - (-a*d + b*c)**7/(3*b**8*(a + b*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1285():
    f = (c + d*x)**7/(a + b*x)**5
    F = 21*d**5*x*(-a*d + b*c)**2/b**7 + d**7*(a + b*x)**3/(3*b**8) + 7*d**6*(a + b*x)**2*(-a*d + b*c)/(2*b**8) + 35*d**4*(-a*d + b*c)**3*log(a + b*x)/b**8 - 35*d**3*(-a*d + b*c)**4/(b**8*(a + b*x)) - 21*d**2*(-a*d + b*c)**5/(2*b**8*(a + b*x)**2) - 7*d*(-a*d + b*c)**6/(3*b**8*(a + b*x)**3) - (-a*d + b*c)**7/(4*b**8*(a + b*x)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1286():
    f = (c + d*x)**7/(a + b*x)**6
    F = d**7*x**2/(2*b**6) + d**6*x*(-6*a*d + 7*b*c)/b**7 + 21*d**5*(-a*d + b*c)**2*log(a + b*x)/b**8 - 35*d**4*(-a*d + b*c)**3/(b**8*(a + b*x)) - 35*d**3*(-a*d + b*c)**4/(2*b**8*(a + b*x)**2) - 7*d**2*(-a*d + b*c)**5/(b**8*(a + b*x)**3) - 7*d*(-a*d + b*c)**6/(4*b**8*(a + b*x)**4) - (-a*d + b*c)**7/(5*b**8*(a + b*x)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1287():
    f = (c + d*x)**7/(a + b*x)**7
    F = d**7*x/b**7 + 7*d**6*(-a*d + b*c)*log(a + b*x)/b**8 - 21*d**5*(-a*d + b*c)**2/(b**8*(a + b*x)) - 35*d**4*(-a*d + b*c)**3/(2*b**8*(a + b*x)**2) - 35*d**3*(-a*d + b*c)**4/(3*b**8*(a + b*x)**3) - 21*d**2*(-a*d + b*c)**5/(4*b**8*(a + b*x)**4) - 7*d*(-a*d + b*c)**6/(5*b**8*(a + b*x)**5) - (-a*d + b*c)**7/(6*b**8*(a + b*x)**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1288():
    f = (c + d*x)**7/(a + b*x)**8
    F = d**7*log(a + b*x)/b**8 - 7*d**6*(-a*d + b*c)/(b**8*(a + b*x)) - 21*d**5*(-a*d + b*c)**2/(2*b**8*(a + b*x)**2) - 35*d**4*(-a*d + b*c)**3/(3*b**8*(a + b*x)**3) - 35*d**3*(-a*d + b*c)**4/(4*b**8*(a + b*x)**4) - 21*d**2*(-a*d + b*c)**5/(5*b**8*(a + b*x)**5) - 7*d*(-a*d + b*c)**6/(6*b**8*(a + b*x)**6) - (-a*d + b*c)**7/(7*b**8*(a + b*x)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1289():
    f = (c + d*x)**7/(a + b*x)**9
    F = -(c + d*x)**8/((a + b*x)**8*(-8*a*d + 8*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1290():
    f = (c + d*x)**7/(a + b*x)**10
    F = d*(c + d*x)**8/(72*(a + b*x)**8*(-a*d + b*c)**2) - (c + d*x)**8/((a + b*x)**9*(-9*a*d + 9*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1291():
    f = (c + d*x)**7/(a + b*x)**11
    F = -d**2*(c + d*x)**8/(360*(a + b*x)**8*(-a*d + b*c)**3) + d*(c + d*x)**8/(45*(a + b*x)**9*(-a*d + b*c)**2) - (c + d*x)**8/((a + b*x)**10*(-10*a*d + 10*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1292():
    f = (c + d*x)**7/(a + b*x)**12
    F = d**3*(c + d*x)**8/(1320*(a + b*x)**8*(-a*d + b*c)**4) - d**2*(c + d*x)**8/(165*(a + b*x)**9*(-a*d + b*c)**3) + 3*d*(c + d*x)**8/(110*(a + b*x)**10*(-a*d + b*c)**2) - (c + d*x)**8/((a + b*x)**11*(-11*a*d + 11*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1293():
    f = (c + d*x)**7/(a + b*x)**13
    F = -d**4*(c + d*x)**8/(3960*(a + b*x)**8*(-a*d + b*c)**5) + d**3*(c + d*x)**8/(495*(a + b*x)**9*(-a*d + b*c)**4) - d**2*(c + d*x)**8/(110*(a + b*x)**10*(-a*d + b*c)**3) + d*(c + d*x)**8/(33*(a + b*x)**11*(-a*d + b*c)**2) - (c + d*x)**8/((a + b*x)**12*(-12*a*d + 12*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1294():
    f = (c + d*x)**7/(a + b*x)**14
    F = -d**7/(6*b**8*(a + b*x)**6) - d**6*(-a*d + b*c)/(b**8*(a + b*x)**7) - 21*d**5*(-a*d + b*c)**2/(8*b**8*(a + b*x)**8) - 35*d**4*(-a*d + b*c)**3/(9*b**8*(a + b*x)**9) - 7*d**3*(-a*d + b*c)**4/(2*b**8*(a + b*x)**10) - 21*d**2*(-a*d + b*c)**5/(11*b**8*(a + b*x)**11) - 7*d*(-a*d + b*c)**6/(12*b**8*(a + b*x)**12) - (-a*d + b*c)**7/(13*b**8*(a + b*x)**13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1295():
    f = (c + d*x)**7/(a + b*x)**15
    F = -d**7/(7*b**8*(a + b*x)**7) - 7*d**6*(-a*d + b*c)/(8*b**8*(a + b*x)**8) - 7*d**5*(-a*d + b*c)**2/(3*b**8*(a + b*x)**9) - 7*d**4*(-a*d + b*c)**3/(2*b**8*(a + b*x)**10) - 35*d**3*(-a*d + b*c)**4/(11*b**8*(a + b*x)**11) - 7*d**2*(-a*d + b*c)**5/(4*b**8*(a + b*x)**12) - 7*d*(-a*d + b*c)**6/(13*b**8*(a + b*x)**13) - (-a*d + b*c)**7/(14*b**8*(a + b*x)**14)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1296():
    f = (c + d*x)**7/(a + b*x)**16
    F = -d**7/(8*b**8*(a + b*x)**8) - 7*d**6*(-a*d + b*c)/(9*b**8*(a + b*x)**9) - 21*d**5*(-a*d + b*c)**2/(10*b**8*(a + b*x)**10) - 35*d**4*(-a*d + b*c)**3/(11*b**8*(a + b*x)**11) - 35*d**3*(-a*d + b*c)**4/(12*b**8*(a + b*x)**12) - 21*d**2*(-a*d + b*c)**5/(13*b**8*(a + b*x)**13) - d*(-a*d + b*c)**6/(2*b**8*(a + b*x)**14) - (-a*d + b*c)**7/(15*b**8*(a + b*x)**15)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1297():
    f = (a + b*x)**12*(c + d*x)**10
    F = d**10*(a + b*x)**23/(23*b**11) + 5*d**9*(a + b*x)**22*(-a*d + b*c)/(11*b**11) + 15*d**8*(a + b*x)**21*(-a*d + b*c)**2/(7*b**11) + 6*d**7*(a + b*x)**20*(-a*d + b*c)**3/b**11 + 210*d**6*(a + b*x)**19*(-a*d + b*c)**4/(19*b**11) + 14*d**5*(a + b*x)**18*(-a*d + b*c)**5/b**11 + 210*d**4*(a + b*x)**17*(-a*d + b*c)**6/(17*b**11) + 15*d**3*(a + b*x)**16*(-a*d + b*c)**7/(2*b**11) + 3*d**2*(a + b*x)**15*(-a*d + b*c)**8/b**11 + 5*d*(a + b*x)**14*(-a*d + b*c)**9/(7*b**11) + (a + b*x)**13*(-a*d + b*c)**10/(13*b**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1298():
    f = (a + b*x)**11*(c + d*x)**10
    F = d**10*(a + b*x)**22/(22*b**11) + 10*d**9*(a + b*x)**21*(-a*d + b*c)/(21*b**11) + 9*d**8*(a + b*x)**20*(-a*d + b*c)**2/(4*b**11) + 120*d**7*(a + b*x)**19*(-a*d + b*c)**3/(19*b**11) + 35*d**6*(a + b*x)**18*(-a*d + b*c)**4/(3*b**11) + 252*d**5*(a + b*x)**17*(-a*d + b*c)**5/(17*b**11) + 105*d**4*(a + b*x)**16*(-a*d + b*c)**6/(8*b**11) + 8*d**3*(a + b*x)**15*(-a*d + b*c)**7/b**11 + 45*d**2*(a + b*x)**14*(-a*d + b*c)**8/(14*b**11) + 10*d*(a + b*x)**13*(-a*d + b*c)**9/(13*b**11) + (a + b*x)**12*(-a*d + b*c)**10/(12*b**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1299():
    f = (a + b*x)**10*(c + d*x)**10
    F = d**10*(a + b*x)**21/(21*b**11) + d**9*(a + b*x)**20*(-a*d + b*c)/(2*b**11) + 45*d**8*(a + b*x)**19*(-a*d + b*c)**2/(19*b**11) + 20*d**7*(a + b*x)**18*(-a*d + b*c)**3/(3*b**11) + 210*d**6*(a + b*x)**17*(-a*d + b*c)**4/(17*b**11) + 63*d**5*(a + b*x)**16*(-a*d + b*c)**5/(4*b**11) + 14*d**4*(a + b*x)**15*(-a*d + b*c)**6/b**11 + 60*d**3*(a + b*x)**14*(-a*d + b*c)**7/(7*b**11) + 45*d**2*(a + b*x)**13*(-a*d + b*c)**8/(13*b**11) + 5*d*(a + b*x)**12*(-a*d + b*c)**9/(6*b**11) + (a + b*x)**11*(-a*d + b*c)**10/(11*b**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1300():
    f = (a + b*x)**9*(c + d*x)**10
    F = b**9*(c + d*x)**20/(20*d**10) - 9*b**8*(c + d*x)**19*(-a*d + b*c)/(19*d**10) + 2*b**7*(c + d*x)**18*(-a*d + b*c)**2/d**10 - 84*b**6*(c + d*x)**17*(-a*d + b*c)**3/(17*d**10) + 63*b**5*(c + d*x)**16*(-a*d + b*c)**4/(8*d**10) - 42*b**4*(c + d*x)**15*(-a*d + b*c)**5/(5*d**10) + 6*b**3*(c + d*x)**14*(-a*d + b*c)**6/d**10 - 36*b**2*(c + d*x)**13*(-a*d + b*c)**7/(13*d**10) + 3*b*(c + d*x)**12*(-a*d + b*c)**8/(4*d**10) - (c + d*x)**11*(-a*d + b*c)**9/(11*d**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1301():
    f = (a + b*x)**8*(c + d*x)**10
    F = b**8*(c + d*x)**19/(19*d**9) - 4*b**7*(c + d*x)**18*(-a*d + b*c)/(9*d**9) + 28*b**6*(c + d*x)**17*(-a*d + b*c)**2/(17*d**9) - 7*b**5*(c + d*x)**16*(-a*d + b*c)**3/(2*d**9) + 14*b**4*(c + d*x)**15*(-a*d + b*c)**4/(3*d**9) - 4*b**3*(c + d*x)**14*(-a*d + b*c)**5/d**9 + 28*b**2*(c + d*x)**13*(-a*d + b*c)**6/(13*d**9) - 2*b*(c + d*x)**12*(-a*d + b*c)**7/(3*d**9) + (c + d*x)**11*(-a*d + b*c)**8/(11*d**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1302():
    f = (a + b*x)**7*(c + d*x)**10
    F = b**7*(c + d*x)**18/(18*d**8) - 7*b**6*(c + d*x)**17*(-a*d + b*c)/(17*d**8) + 21*b**5*(c + d*x)**16*(-a*d + b*c)**2/(16*d**8) - 7*b**4*(c + d*x)**15*(-a*d + b*c)**3/(3*d**8) + 5*b**3*(c + d*x)**14*(-a*d + b*c)**4/(2*d**8) - 21*b**2*(c + d*x)**13*(-a*d + b*c)**5/(13*d**8) + 7*b*(c + d*x)**12*(-a*d + b*c)**6/(12*d**8) - (c + d*x)**11*(-a*d + b*c)**7/(11*d**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1303():
    f = (a + b*x)**6*(c + d*x)**10
    F = b**6*(c + d*x)**17/(17*d**7) - 3*b**5*(c + d*x)**16*(-a*d + b*c)/(8*d**7) + b**4*(c + d*x)**15*(-a*d + b*c)**2/d**7 - 10*b**3*(c + d*x)**14*(-a*d + b*c)**3/(7*d**7) + 15*b**2*(c + d*x)**13*(-a*d + b*c)**4/(13*d**7) - b*(c + d*x)**12*(-a*d + b*c)**5/(2*d**7) + (c + d*x)**11*(-a*d + b*c)**6/(11*d**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1304():
    f = (a + b*x)**5*(c + d*x)**10
    F = b**5*(c + d*x)**16/(16*d**6) - b**4*(c + d*x)**15*(-a*d + b*c)/(3*d**6) + 5*b**3*(c + d*x)**14*(-a*d + b*c)**2/(7*d**6) - 10*b**2*(c + d*x)**13*(-a*d + b*c)**3/(13*d**6) + 5*b*(c + d*x)**12*(-a*d + b*c)**4/(12*d**6) - (c + d*x)**11*(-a*d + b*c)**5/(11*d**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1305():
    f = (a + b*x)**4*(c + d*x)**10
    F = b**4*(c + d*x)**15/(15*d**5) - 2*b**3*(c + d*x)**14*(-a*d + b*c)/(7*d**5) + 6*b**2*(c + d*x)**13*(-a*d + b*c)**2/(13*d**5) - b*(c + d*x)**12*(-a*d + b*c)**3/(3*d**5) + (c + d*x)**11*(-a*d + b*c)**4/(11*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1306():
    f = (a + b*x)**3*(c + d*x)**10
    F = b**3*(c + d*x)**14/(14*d**4) - 3*b**2*(c + d*x)**13*(-a*d + b*c)/(13*d**4) + b*(c + d*x)**12*(-a*d + b*c)**2/(4*d**4) - (c + d*x)**11*(-a*d + b*c)**3/(11*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1307():
    f = (a + b*x)**2*(c + d*x)**10
    F = b**2*(c + d*x)**13/(13*d**3) - b*(c + d*x)**12*(-a*d + b*c)/(6*d**3) + (c + d*x)**11*(-a*d + b*c)**2/(11*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1308():
    f = (a + b*x)*(c + d*x)**10
    F = b*(c + d*x)**12/(12*d**2) - (c + d*x)**11*(-a*d + b*c)/(11*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1309():
    f = (c + d*x)**10
    F = (c + d*x)**11/(11*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1310():
    f = (c + d*x)**10/(a + b*x)
    F = (c + d*x)**10/(10*b) + (c + d*x)**9*(-a*d + b*c)/(9*b**2) + (c + d*x)**8*(-a*d + b*c)**2/(8*b**3) + (c + d*x)**7*(-a*d + b*c)**3/(7*b**4) + (c + d*x)**6*(-a*d + b*c)**4/(6*b**5) + (c + d*x)**5*(-a*d + b*c)**5/(5*b**6) + (c + d*x)**4*(-a*d + b*c)**6/(4*b**7) + (c + d*x)**3*(-a*d + b*c)**7/(3*b**8) + (c + d*x)**2*(-a*d + b*c)**8/(2*b**9) + d*x*(-a*d + b*c)**9/b**10 + (-a*d + b*c)**10*log(a + b*x)/b**11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1311():
    f = (c + d*x)**10/(a + b*x)**2
    F = 45*d**2*x*(-a*d + b*c)**8/b**10 + d**10*(a + b*x)**9/(9*b**11) + 5*d**9*(a + b*x)**8*(-a*d + b*c)/(4*b**11) + 45*d**8*(a + b*x)**7*(-a*d + b*c)**2/(7*b**11) + 20*d**7*(a + b*x)**6*(-a*d + b*c)**3/b**11 + 42*d**6*(a + b*x)**5*(-a*d + b*c)**4/b**11 + 63*d**5*(a + b*x)**4*(-a*d + b*c)**5/b**11 + 70*d**4*(a + b*x)**3*(-a*d + b*c)**6/b**11 + 60*d**3*(a + b*x)**2*(-a*d + b*c)**7/b**11 + 10*d*(-a*d + b*c)**9*log(a + b*x)/b**11 - (-a*d + b*c)**10/(b**11*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1312():
    f = (c + d*x)**10/(a + b*x)**3
    F = 120*d**3*x*(-a*d + b*c)**7/b**10 + d**10*(a + b*x)**8/(8*b**11) + 10*d**9*(a + b*x)**7*(-a*d + b*c)/(7*b**11) + 15*d**8*(a + b*x)**6*(-a*d + b*c)**2/(2*b**11) + 24*d**7*(a + b*x)**5*(-a*d + b*c)**3/b**11 + 105*d**6*(a + b*x)**4*(-a*d + b*c)**4/(2*b**11) + 84*d**5*(a + b*x)**3*(-a*d + b*c)**5/b**11 + 105*d**4*(a + b*x)**2*(-a*d + b*c)**6/b**11 + 45*d**2*(-a*d + b*c)**8*log(a + b*x)/b**11 - 10*d*(-a*d + b*c)**9/(b**11*(a + b*x)) - (-a*d + b*c)**10/(2*b**11*(a + b*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1313():
    f = (c + d*x)**10/(a + b*x)**4
    F = 210*d**4*x*(-a*d + b*c)**6/b**10 + d**10*(a + b*x)**7/(7*b**11) + 5*d**9*(a + b*x)**6*(-a*d + b*c)/(3*b**11) + 9*d**8*(a + b*x)**5*(-a*d + b*c)**2/b**11 + 30*d**7*(a + b*x)**4*(-a*d + b*c)**3/b**11 + 70*d**6*(a + b*x)**3*(-a*d + b*c)**4/b**11 + 126*d**5*(a + b*x)**2*(-a*d + b*c)**5/b**11 + 120*d**3*(-a*d + b*c)**7*log(a + b*x)/b**11 - 45*d**2*(-a*d + b*c)**8/(b**11*(a + b*x)) - 5*d*(-a*d + b*c)**9/(b**11*(a + b*x)**2) - (-a*d + b*c)**10/(3*b**11*(a + b*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1314():
    f = (c + d*x)**10/(a + b*x)**5
    F = 252*d**5*x*(-a*d + b*c)**5/b**10 + d**10*(a + b*x)**6/(6*b**11) + 2*d**9*(a + b*x)**5*(-a*d + b*c)/b**11 + 45*d**8*(a + b*x)**4*(-a*d + b*c)**2/(4*b**11) + 40*d**7*(a + b*x)**3*(-a*d + b*c)**3/b**11 + 105*d**6*(a + b*x)**2*(-a*d + b*c)**4/b**11 + 210*d**4*(-a*d + b*c)**6*log(a + b*x)/b**11 - 120*d**3*(-a*d + b*c)**7/(b**11*(a + b*x)) - 45*d**2*(-a*d + b*c)**8/(2*b**11*(a + b*x)**2) - 10*d*(-a*d + b*c)**9/(3*b**11*(a + b*x)**3) - (-a*d + b*c)**10/(4*b**11*(a + b*x)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1315():
    f = (c + d*x)**10/(a + b*x)**6
    F = 210*d**6*x*(-a*d + b*c)**4/b**10 + d**10*(a + b*x)**5/(5*b**11) + 5*d**9*(a + b*x)**4*(-a*d + b*c)/(2*b**11) + 15*d**8*(a + b*x)**3*(-a*d + b*c)**2/b**11 + 60*d**7*(a + b*x)**2*(-a*d + b*c)**3/b**11 + 252*d**5*(-a*d + b*c)**5*log(a + b*x)/b**11 - 210*d**4*(-a*d + b*c)**6/(b**11*(a + b*x)) - 60*d**3*(-a*d + b*c)**7/(b**11*(a + b*x)**2) - 15*d**2*(-a*d + b*c)**8/(b**11*(a + b*x)**3) - 5*d*(-a*d + b*c)**9/(2*b**11*(a + b*x)**4) - (-a*d + b*c)**10/(5*b**11*(a + b*x)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1316():
    f = (c + d*x)**10/(a + b*x)**7
    F = 120*d**7*x*(-a*d + b*c)**3/b**10 + d**10*(a + b*x)**4/(4*b**11) + 10*d**9*(a + b*x)**3*(-a*d + b*c)/(3*b**11) + 45*d**8*(a + b*x)**2*(-a*d + b*c)**2/(2*b**11) + 210*d**6*(-a*d + b*c)**4*log(a + b*x)/b**11 - 252*d**5*(-a*d + b*c)**5/(b**11*(a + b*x)) - 105*d**4*(-a*d + b*c)**6/(b**11*(a + b*x)**2) - 40*d**3*(-a*d + b*c)**7/(b**11*(a + b*x)**3) - 45*d**2*(-a*d + b*c)**8/(4*b**11*(a + b*x)**4) - 2*d*(-a*d + b*c)**9/(b**11*(a + b*x)**5) - (-a*d + b*c)**10/(6*b**11*(a + b*x)**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1317():
    f = (c + d*x)**10/(a + b*x)**8
    F = 45*d**8*x*(-a*d + b*c)**2/b**10 + d**10*(a + b*x)**3/(3*b**11) + 5*d**9*(a + b*x)**2*(-a*d + b*c)/b**11 + 120*d**7*(-a*d + b*c)**3*log(a + b*x)/b**11 - 210*d**6*(-a*d + b*c)**4/(b**11*(a + b*x)) - 126*d**5*(-a*d + b*c)**5/(b**11*(a + b*x)**2) - 70*d**4*(-a*d + b*c)**6/(b**11*(a + b*x)**3) - 30*d**3*(-a*d + b*c)**7/(b**11*(a + b*x)**4) - 9*d**2*(-a*d + b*c)**8/(b**11*(a + b*x)**5) - 5*d*(-a*d + b*c)**9/(3*b**11*(a + b*x)**6) - (-a*d + b*c)**10/(7*b**11*(a + b*x)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1318():
    f = (c + d*x)**10/(a + b*x)**9
    F = d**10*x**2/(2*b**9) + d**9*x*(-9*a*d + 10*b*c)/b**10 + 45*d**8*(-a*d + b*c)**2*log(a + b*x)/b**11 - 120*d**7*(-a*d + b*c)**3/(b**11*(a + b*x)) - 105*d**6*(-a*d + b*c)**4/(b**11*(a + b*x)**2) - 84*d**5*(-a*d + b*c)**5/(b**11*(a + b*x)**3) - 105*d**4*(-a*d + b*c)**6/(2*b**11*(a + b*x)**4) - 24*d**3*(-a*d + b*c)**7/(b**11*(a + b*x)**5) - 15*d**2*(-a*d + b*c)**8/(2*b**11*(a + b*x)**6) - 10*d*(-a*d + b*c)**9/(7*b**11*(a + b*x)**7) - (-a*d + b*c)**10/(8*b**11*(a + b*x)**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1319():
    f = (c + d*x)**10/(a + b*x)**10
    F = d**10*x/b**10 + 10*d**9*(-a*d + b*c)*log(a + b*x)/b**11 - 45*d**8*(-a*d + b*c)**2/(b**11*(a + b*x)) - 60*d**7*(-a*d + b*c)**3/(b**11*(a + b*x)**2) - 70*d**6*(-a*d + b*c)**4/(b**11*(a + b*x)**3) - 63*d**5*(-a*d + b*c)**5/(b**11*(a + b*x)**4) - 42*d**4*(-a*d + b*c)**6/(b**11*(a + b*x)**5) - 20*d**3*(-a*d + b*c)**7/(b**11*(a + b*x)**6) - 45*d**2*(-a*d + b*c)**8/(7*b**11*(a + b*x)**7) - 5*d*(-a*d + b*c)**9/(4*b**11*(a + b*x)**8) - (-a*d + b*c)**10/(9*b**11*(a + b*x)**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1320():
    f = (c + d*x)**10/(a + b*x)**11
    F = d**10*log(a + b*x)/b**11 - 10*d**9*(-a*d + b*c)/(b**11*(a + b*x)) - 45*d**8*(-a*d + b*c)**2/(2*b**11*(a + b*x)**2) - 40*d**7*(-a*d + b*c)**3/(b**11*(a + b*x)**3) - 105*d**6*(-a*d + b*c)**4/(2*b**11*(a + b*x)**4) - 252*d**5*(-a*d + b*c)**5/(5*b**11*(a + b*x)**5) - 35*d**4*(-a*d + b*c)**6/(b**11*(a + b*x)**6) - 120*d**3*(-a*d + b*c)**7/(7*b**11*(a + b*x)**7) - 45*d**2*(-a*d + b*c)**8/(8*b**11*(a + b*x)**8) - 10*d*(-a*d + b*c)**9/(9*b**11*(a + b*x)**9) - (-a*d + b*c)**10/(10*b**11*(a + b*x)**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1321():
    f = (c + d*x)**10/(a + b*x)**12
    F = -(c + d*x)**11/((a + b*x)**11*(-11*a*d + 11*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1322():
    f = (c + d*x)**10/(a + b*x)**13
    F = d*(c + d*x)**11/(132*(a + b*x)**11*(-a*d + b*c)**2) - (c + d*x)**11/((a + b*x)**12*(-12*a*d + 12*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1323():
    f = (c + d*x)**10/(a + b*x)**14
    F = -d**2*(c + d*x)**11/(858*(a + b*x)**11*(-a*d + b*c)**3) + d*(c + d*x)**11/(78*(a + b*x)**12*(-a*d + b*c)**2) - (c + d*x)**11/((a + b*x)**13*(-13*a*d + 13*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1324():
    f = (c + d*x)**10/(a + b*x)**15
    F = d**3*(c + d*x)**11/(4004*(a + b*x)**11*(-a*d + b*c)**4) - d**2*(c + d*x)**11/(364*(a + b*x)**12*(-a*d + b*c)**3) + 3*d*(c + d*x)**11/(182*(a + b*x)**13*(-a*d + b*c)**2) - (c + d*x)**11/((a + b*x)**14*(-14*a*d + 14*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1325():
    f = (c + d*x)**10/(a + b*x)**16
    F = -d**4*(c + d*x)**11/(15015*(a + b*x)**11*(-a*d + b*c)**5) + d**3*(c + d*x)**11/(1365*(a + b*x)**12*(-a*d + b*c)**4) - 2*d**2*(c + d*x)**11/(455*(a + b*x)**13*(-a*d + b*c)**3) + 2*d*(c + d*x)**11/(105*(a + b*x)**14*(-a*d + b*c)**2) - (c + d*x)**11/((a + b*x)**15*(-15*a*d + 15*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1326():
    f = (c + d*x)**10/(a + b*x)**17
    F = d**5*(c + d*x)**11/(48048*(a + b*x)**11*(-a*d + b*c)**6) - d**4*(c + d*x)**11/(4368*(a + b*x)**12*(-a*d + b*c)**5) + d**3*(c + d*x)**11/(728*(a + b*x)**13*(-a*d + b*c)**4) - d**2*(c + d*x)**11/(168*(a + b*x)**14*(-a*d + b*c)**3) + d*(c + d*x)**11/(48*(a + b*x)**15*(-a*d + b*c)**2) - (c + d*x)**11/((a + b*x)**16*(-16*a*d + 16*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1327():
    f = (c + d*x)**10/(a + b*x)**18
    F = -d**6*(c + d*x)**11/(136136*(a + b*x)**11*(-a*d + b*c)**7) + d**5*(c + d*x)**11/(12376*(a + b*x)**12*(-a*d + b*c)**6) - 3*d**4*(c + d*x)**11/(6188*(a + b*x)**13*(-a*d + b*c)**5) + d**3*(c + d*x)**11/(476*(a + b*x)**14*(-a*d + b*c)**4) - d**2*(c + d*x)**11/(136*(a + b*x)**15*(-a*d + b*c)**3) + 3*d*(c + d*x)**11/(136*(a + b*x)**16*(-a*d + b*c)**2) - (c + d*x)**11/((a + b*x)**17*(-17*a*d + 17*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1328():
    f = (c + d*x)**10/(a + b*x)**19
    F = d**7*(c + d*x)**11/(350064*(a + b*x)**11*(-a*d + b*c)**8) - d**6*(c + d*x)**11/(31824*(a + b*x)**12*(-a*d + b*c)**7) + d**5*(c + d*x)**11/(5304*(a + b*x)**13*(-a*d + b*c)**6) - d**4*(c + d*x)**11/(1224*(a + b*x)**14*(-a*d + b*c)**5) + 7*d**3*(c + d*x)**11/(2448*(a + b*x)**15*(-a*d + b*c)**4) - 7*d**2*(c + d*x)**11/(816*(a + b*x)**16*(-a*d + b*c)**3) + 7*d*(c + d*x)**11/(306*(a + b*x)**17*(-a*d + b*c)**2) - (c + d*x)**11/((a + b*x)**18*(-18*a*d + 18*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1329():
    f = (c + d*x)**10/(a + b*x)**20
    F = -d**10/(9*b**11*(a + b*x)**9) - d**9*(-a*d + b*c)/(b**11*(a + b*x)**10) - 45*d**8*(-a*d + b*c)**2/(11*b**11*(a + b*x)**11) - 10*d**7*(-a*d + b*c)**3/(b**11*(a + b*x)**12) - 210*d**6*(-a*d + b*c)**4/(13*b**11*(a + b*x)**13) - 18*d**5*(-a*d + b*c)**5/(b**11*(a + b*x)**14) - 14*d**4*(-a*d + b*c)**6/(b**11*(a + b*x)**15) - 15*d**3*(-a*d + b*c)**7/(2*b**11*(a + b*x)**16) - 45*d**2*(-a*d + b*c)**8/(17*b**11*(a + b*x)**17) - 5*d*(-a*d + b*c)**9/(9*b**11*(a + b*x)**18) - (-a*d + b*c)**10/(19*b**11*(a + b*x)**19)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1330():
    f = (c + d*x)**10/(a + b*x)**21
    F = -d**10/(10*b**11*(a + b*x)**10) - 10*d**9*(-a*d + b*c)/(11*b**11*(a + b*x)**11) - 15*d**8*(-a*d + b*c)**2/(4*b**11*(a + b*x)**12) - 120*d**7*(-a*d + b*c)**3/(13*b**11*(a + b*x)**13) - 15*d**6*(-a*d + b*c)**4/(b**11*(a + b*x)**14) - 84*d**5*(-a*d + b*c)**5/(5*b**11*(a + b*x)**15) - 105*d**4*(-a*d + b*c)**6/(8*b**11*(a + b*x)**16) - 120*d**3*(-a*d + b*c)**7/(17*b**11*(a + b*x)**17) - 5*d**2*(-a*d + b*c)**8/(2*b**11*(a + b*x)**18) - 10*d*(-a*d + b*c)**9/(19*b**11*(a + b*x)**19) - (-a*d + b*c)**10/(20*b**11*(a + b*x)**20)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1331():
    f = (c + d*x)**10/(a + b*x)**22
    F = -d**10/(11*b**11*(a + b*x)**11) - 5*d**9*(-a*d + b*c)/(6*b**11*(a + b*x)**12) - 45*d**8*(-a*d + b*c)**2/(13*b**11*(a + b*x)**13) - 60*d**7*(-a*d + b*c)**3/(7*b**11*(a + b*x)**14) - 14*d**6*(-a*d + b*c)**4/(b**11*(a + b*x)**15) - 63*d**5*(-a*d + b*c)**5/(4*b**11*(a + b*x)**16) - 210*d**4*(-a*d + b*c)**6/(17*b**11*(a + b*x)**17) - 20*d**3*(-a*d + b*c)**7/(3*b**11*(a + b*x)**18) - 45*d**2*(-a*d + b*c)**8/(19*b**11*(a + b*x)**19) - d*(-a*d + b*c)**9/(2*b**11*(a + b*x)**20) - (-a*d + b*c)**10/(21*b**11*(a + b*x)**21)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1332():
    f = (a + b*x)**5/(c + d*x)
    F = b*x*(-a*d + b*c)**4/d**5 + (a + b*x)**5/(5*d) - (a + b*x)**4*(-a*d + b*c)/(4*d**2) + (a + b*x)**3*(-a*d + b*c)**2/(3*d**3) - (a + b*x)**2*(-a*d + b*c)**3/(2*d**4) - (-a*d + b*c)**5*log(c + d*x)/d**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1333():
    f = (a + b*x)**4/(c + d*x)
    F = -b*x*(-a*d + b*c)**3/d**4 + (a + b*x)**4/(4*d) - (a + b*x)**3*(-a*d + b*c)/(3*d**2) + (a + b*x)**2*(-a*d + b*c)**2/(2*d**3) + (-a*d + b*c)**4*log(c + d*x)/d**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1334():
    f = (a + b*x)**3/(c + d*x)
    F = b*x*(-a*d + b*c)**2/d**3 + (a + b*x)**3/(3*d) - (a + b*x)**2*(-a*d + b*c)/(2*d**2) - (-a*d + b*c)**3*log(c + d*x)/d**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1335():
    f = (a + b*x)**2/(c + d*x)
    F = -b*x*(-a*d + b*c)/d**2 + (a + b*x)**2/(2*d) + (-a*d + b*c)**2*log(c + d*x)/d**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1336():
    f = (a + b*x)/(c + d*x)
    F = b*x/d - (-a*d + b*c)*log(c + d*x)/d**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1337():
    f = 1/(c + d*x)
    F = log(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1338():
    f = 1/((a + b*x)*(c + d*x))
    F = log(a + b*x)/(-a*d + b*c) - log(c + d*x)/(-a*d + b*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1339():
    f = 1/((a + b*x)**2*(c + d*x))
    F = -d*log(a + b*x)/(-a*d + b*c)**2 + d*log(c + d*x)/(-a*d + b*c)**2 - 1/((a + b*x)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1340():
    f = 1/((a + b*x)**3*(c + d*x))
    F = d**2*log(a + b*x)/(-a*d + b*c)**3 - d**2*log(c + d*x)/(-a*d + b*c)**3 + d/((a + b*x)*(-a*d + b*c)**2) - 1/((a + b*x)**2*(-2*a*d + 2*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1341():
    f = (a + b*x)**5/(c + d*x)**2
    F = b**5*(c + d*x)**4/(4*d**6) - 5*b**4*(c + d*x)**3*(-a*d + b*c)/(3*d**6) + 5*b**3*(c + d*x)**2*(-a*d + b*c)**2/d**6 - 10*b**2*x*(-a*d + b*c)**3/d**5 + 5*b*(-a*d + b*c)**4*log(c + d*x)/d**6 + (-a*d + b*c)**5/(d**6*(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1342():
    f = (a + b*x)**4/(c + d*x)**2
    F = b**4*(c + d*x)**3/(3*d**5) - 2*b**3*(c + d*x)**2*(-a*d + b*c)/d**5 + 6*b**2*x*(-a*d + b*c)**2/d**4 - 4*b*(-a*d + b*c)**3*log(c + d*x)/d**5 - (-a*d + b*c)**4/(d**5*(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1343():
    f = (a + b*x)**3/(c + d*x)**2
    F = b**3*x**2/(2*d**2) - b**2*x*(-3*a*d + 2*b*c)/d**3 + 3*b*(-a*d + b*c)**2*log(c + d*x)/d**4 + (-a*d + b*c)**3/(d**4*(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1344():
    f = (a + b*x)**2/(c + d*x)**2
    F = b**2*x/d**2 - 2*b*(-a*d + b*c)*log(c + d*x)/d**3 - (-a*d + b*c)**2/(d**3*(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1345():
    f = (a + b*x)/(c + d*x)**2
    F = b*log(c + d*x)/d**2 + (-a*d + b*c)/(d**2*(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1346():
    f = (c + d*x)**(-2)
    F = -1/(d*(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1347():
    f = 1/((a + b*x)*(c + d*x)**2)
    F = b*log(a + b*x)/(-a*d + b*c)**2 - b*log(c + d*x)/(-a*d + b*c)**2 + 1/((c + d*x)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1348():
    f = 1/((a + b*x)**2*(c + d*x)**2)
    F = -2*b*d*log(a + b*x)/(-a*d + b*c)**3 + 2*b*d*log(c + d*x)/(-a*d + b*c)**3 - b/((a + b*x)*(-a*d + b*c)**2) - d/((c + d*x)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1349():
    f = 1/((a + b*x)**3*(c + d*x)**2)
    F = 3*b*d**2*log(a + b*x)/(-a*d + b*c)**4 - 3*b*d**2*log(c + d*x)/(-a*d + b*c)**4 + 2*b*d/((a + b*x)*(-a*d + b*c)**3) - b/(2*(a + b*x)**2*(-a*d + b*c)**2) + d**2/((c + d*x)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1350():
    f = (a + b*x)**6/(c + d*x)**3
    F = b**6*(c + d*x)**4/(4*d**7) - 2*b**5*(c + d*x)**3*(-a*d + b*c)/d**7 + 15*b**4*(c + d*x)**2*(-a*d + b*c)**2/(2*d**7) - 20*b**3*x*(-a*d + b*c)**3/d**6 + 15*b**2*(-a*d + b*c)**4*log(c + d*x)/d**7 + 6*b*(-a*d + b*c)**5/(d**7*(c + d*x)) - (-a*d + b*c)**6/(2*d**7*(c + d*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1351():
    f = (a + b*x)**5/(c + d*x)**3
    F = b**5*(c + d*x)**3/(3*d**6) - 5*b**4*(c + d*x)**2*(-a*d + b*c)/(2*d**6) + 10*b**3*x*(-a*d + b*c)**2/d**5 - 10*b**2*(-a*d + b*c)**3*log(c + d*x)/d**6 - 5*b*(-a*d + b*c)**4/(d**6*(c + d*x)) + (-a*d + b*c)**5/(2*d**6*(c + d*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1352():
    f = (a + b*x)**4/(c + d*x)**3
    F = b**4*x**2/(2*d**3) - b**3*x*(-4*a*d + 3*b*c)/d**4 + 6*b**2*(-a*d + b*c)**2*log(c + d*x)/d**5 + 4*b*(-a*d + b*c)**3/(d**5*(c + d*x)) - (-a*d + b*c)**4/(2*d**5*(c + d*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1353():
    f = (a + b*x)**3/(c + d*x)**3
    F = b**3*x/d**3 - 3*b**2*(-a*d + b*c)*log(c + d*x)/d**4 - 3*b*(-a*d + b*c)**2/(d**4*(c + d*x)) + (-a*d + b*c)**3/(2*d**4*(c + d*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1354():
    f = (a + b*x)**2/(c + d*x)**3
    F = b**2*log(c + d*x)/d**3 + 2*b*(-a*d + b*c)/(d**3*(c + d*x)) - (-a*d + b*c)**2/(2*d**3*(c + d*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1355():
    f = (a + b*x)/(c + d*x)**3
    F = (a + b*x)**2/((c + d*x)**2*(-2*a*d + 2*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1356():
    f = (c + d*x)**(-3)
    F = -1/(2*d*(c + d*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1357():
    f = 1/((a + b*x)*(c + d*x)**3)
    F = b**2*log(a + b*x)/(-a*d + b*c)**3 - b**2*log(c + d*x)/(-a*d + b*c)**3 + b/((c + d*x)*(-a*d + b*c)**2) + 1/((c + d*x)**2*(-2*a*d + 2*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1358():
    f = 1/((a + b*x)**2*(c + d*x)**3)
    F = -3*b**2*d*log(a + b*x)/(-a*d + b*c)**4 + 3*b**2*d*log(c + d*x)/(-a*d + b*c)**4 - b**2/((a + b*x)*(-a*d + b*c)**3) - 2*b*d/((c + d*x)*(-a*d + b*c)**3) - d/(2*(c + d*x)**2*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1359():
    f = 1/((a + b*x)**3*(c + d*x)**3)
    F = 6*b**2*d**2*log(a + b*x)/(-a*d + b*c)**5 - 6*b**2*d**2*log(c + d*x)/(-a*d + b*c)**5 + 3*b**2*d/((a + b*x)*(-a*d + b*c)**4) - b**2/(2*(a + b*x)**2*(-a*d + b*c)**3) + 3*b*d**2/((c + d*x)*(-a*d + b*c)**4) + d**2/(2*(c + d*x)**2*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1360():
    f = (a + b*x)**9/(c + d*x)**8
    F = b**9*x**2/(2*d**8) - b**8*x*(-9*a*d + 8*b*c)/d**9 + 36*b**7*(-a*d + b*c)**2*log(c + d*x)/d**10 + 84*b**6*(-a*d + b*c)**3/(d**10*(c + d*x)) - 63*b**5*(-a*d + b*c)**4/(d**10*(c + d*x)**2) + 42*b**4*(-a*d + b*c)**5/(d**10*(c + d*x)**3) - 21*b**3*(-a*d + b*c)**6/(d**10*(c + d*x)**4) + 36*b**2*(-a*d + b*c)**7/(5*d**10*(c + d*x)**5) - 3*b*(-a*d + b*c)**8/(2*d**10*(c + d*x)**6) + (-a*d + b*c)**9/(7*d**10*(c + d*x)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1361():
    f = (a + b*x)**8/(c + d*x)**8
    F = b**8*x/d**8 - 8*b**7*(-a*d + b*c)*log(c + d*x)/d**9 - 28*b**6*(-a*d + b*c)**2/(d**9*(c + d*x)) + 28*b**5*(-a*d + b*c)**3/(d**9*(c + d*x)**2) - 70*b**4*(-a*d + b*c)**4/(3*d**9*(c + d*x)**3) + 14*b**3*(-a*d + b*c)**5/(d**9*(c + d*x)**4) - 28*b**2*(-a*d + b*c)**6/(5*d**9*(c + d*x)**5) + 4*b*(-a*d + b*c)**7/(3*d**9*(c + d*x)**6) - (-a*d + b*c)**8/(7*d**9*(c + d*x)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1362():
    f = (a + b*x)**7/(c + d*x)**8
    F = b**7*log(c + d*x)/d**8 + 7*b**6*(-a*d + b*c)/(d**8*(c + d*x)) - 21*b**5*(-a*d + b*c)**2/(2*d**8*(c + d*x)**2) + 35*b**4*(-a*d + b*c)**3/(3*d**8*(c + d*x)**3) - 35*b**3*(-a*d + b*c)**4/(4*d**8*(c + d*x)**4) + 21*b**2*(-a*d + b*c)**5/(5*d**8*(c + d*x)**5) - 7*b*(-a*d + b*c)**6/(6*d**8*(c + d*x)**6) + (-a*d + b*c)**7/(7*d**8*(c + d*x)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1363():
    f = (a + b*x)**6/(c + d*x)**8
    F = (a + b*x)**7/((c + d*x)**7*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1364():
    f = (a + b*x)**5/(c + d*x)**8
    F = b*(a + b*x)**6/(42*(c + d*x)**6*(-a*d + b*c)**2) + (a + b*x)**6/((c + d*x)**7*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1365():
    f = (a + b*x)**4/(c + d*x)**8
    F = b**2*(a + b*x)**5/(105*(c + d*x)**5*(-a*d + b*c)**3) + b*(a + b*x)**5/(21*(c + d*x)**6*(-a*d + b*c)**2) + (a + b*x)**5/((c + d*x)**7*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1366():
    f = (a + b*x)**3/(c + d*x)**8
    F = -b**3/(4*d**4*(c + d*x)**4) + 3*b**2*(-a*d + b*c)/(5*d**4*(c + d*x)**5) - b*(-a*d + b*c)**2/(2*d**4*(c + d*x)**6) + (-a*d + b*c)**3/(7*d**4*(c + d*x)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1367():
    f = (a + b*x)**2/(c + d*x)**8
    F = -b**2/(5*d**3*(c + d*x)**5) + b*(-a*d + b*c)/(3*d**3*(c + d*x)**6) - (-a*d + b*c)**2/(7*d**3*(c + d*x)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1368():
    f = (a + b*x)/(c + d*x)**8
    F = -b/(6*d**2*(c + d*x)**6) + (-a*d + b*c)/(7*d**2*(c + d*x)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1369():
    f = (c + d*x)**(-8)
    F = -1/(7*d*(c + d*x)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1370():
    f = 1/((a + b*x)*(c + d*x)**8)
    F = b**7*log(a + b*x)/(-a*d + b*c)**8 - b**7*log(c + d*x)/(-a*d + b*c)**8 + b**6/((c + d*x)*(-a*d + b*c)**7) + b**5/(2*(c + d*x)**2*(-a*d + b*c)**6) + b**4/(3*(c + d*x)**3*(-a*d + b*c)**5) + b**3/(4*(c + d*x)**4*(-a*d + b*c)**4) + b**2/(5*(c + d*x)**5*(-a*d + b*c)**3) + b/(6*(c + d*x)**6*(-a*d + b*c)**2) + 1/((c + d*x)**7*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1371():
    f = 1/((a + b*x)**2*(c + d*x)**8)
    F = -8*b**7*d*log(a + b*x)/(-a*d + b*c)**9 + 8*b**7*d*log(c + d*x)/(-a*d + b*c)**9 - b**7/((a + b*x)*(-a*d + b*c)**8) - 7*b**6*d/((c + d*x)*(-a*d + b*c)**8) - 3*b**5*d/((c + d*x)**2*(-a*d + b*c)**7) - 5*b**4*d/(3*(c + d*x)**3*(-a*d + b*c)**6) - b**3*d/((c + d*x)**4*(-a*d + b*c)**5) - 3*b**2*d/(5*(c + d*x)**5*(-a*d + b*c)**4) - b*d/(3*(c + d*x)**6*(-a*d + b*c)**3) - d/(7*(c + d*x)**7*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1372():
    f = 1/((a + b*x)**3*(c + d*x)**8)
    F = 36*b**7*d**2*log(a + b*x)/(-a*d + b*c)**10 - 36*b**7*d**2*log(c + d*x)/(-a*d + b*c)**10 + 8*b**7*d/((a + b*x)*(-a*d + b*c)**9) - b**7/(2*(a + b*x)**2*(-a*d + b*c)**8) + 28*b**6*d**2/((c + d*x)*(-a*d + b*c)**9) + 21*b**5*d**2/(2*(c + d*x)**2*(-a*d + b*c)**8) + 5*b**4*d**2/((c + d*x)**3*(-a*d + b*c)**7) + 5*b**3*d**2/(2*(c + d*x)**4*(-a*d + b*c)**6) + 6*b**2*d**2/(5*(c + d*x)**5*(-a*d + b*c)**5) + b*d**2/(2*(c + d*x)**6*(-a*d + b*c)**4) + d**2/(7*(c + d*x)**7*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1373():
    f = (a + b*x)**5*sqrt(c + d*x)
    F = 2*b**5*(c + d*x)**(sympy.S(13)/2)/(13*d**6) - 10*b**4*(c + d*x)**(sympy.S(11)/2)*(-a*d + b*c)/(11*d**6) + 20*b**3*(c + d*x)**(sympy.S(9)/2)*(-a*d + b*c)**2/(9*d**6) - 20*b**2*(c + d*x)**(sympy.S(7)/2)*(-a*d + b*c)**3/(7*d**6) + 2*b*(c + d*x)**(sympy.S(5)/2)*(-a*d + b*c)**4/d**6 - 2*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**5/(3*d**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1374():
    f = (a + b*x)**4*sqrt(c + d*x)
    F = 2*b**4*(c + d*x)**(sympy.S(11)/2)/(11*d**5) - 8*b**3*(c + d*x)**(sympy.S(9)/2)*(-a*d + b*c)/(9*d**5) + 12*b**2*(c + d*x)**(sympy.S(7)/2)*(-a*d + b*c)**2/(7*d**5) - 8*b*(c + d*x)**(sympy.S(5)/2)*(-a*d + b*c)**3/(5*d**5) + 2*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**4/(3*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1375():
    f = (a + b*x)**3*sqrt(c + d*x)
    F = 2*b**3*(c + d*x)**(sympy.S(9)/2)/(9*d**4) - 6*b**2*(c + d*x)**(sympy.S(7)/2)*(-a*d + b*c)/(7*d**4) + 6*b*(c + d*x)**(sympy.S(5)/2)*(-a*d + b*c)**2/(5*d**4) - 2*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**3/(3*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1376():
    f = (a + b*x)**2*sqrt(c + d*x)
    F = 2*b**2*(c + d*x)**(sympy.S(7)/2)/(7*d**3) - 4*b*(c + d*x)**(sympy.S(5)/2)*(-a*d + b*c)/(5*d**3) + 2*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**2/(3*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1377():
    f = (a + b*x)*sqrt(c + d*x)
    F = 2*b*(c + d*x)**(sympy.S(5)/2)/(5*d**2) - (c + d*x)**(sympy.S(3)/2)*(-2*a*d + 2*b*c)/(3*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1378():
    f = sqrt(c + d*x)
    F = 2*(c + d*x)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1379():
    f = sqrt(c + d*x)/(a + b*x)
    F = 2*sqrt(c + d*x)/b - 2*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1380():
    f = sqrt(c + d*x)/(a + b*x)**2
    F = -sqrt(c + d*x)/(b*(a + b*x)) - d*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(b**(sympy.S(3)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1381():
    f = sqrt(c + d*x)/(a + b*x)**3
    F = -d*sqrt(c + d*x)/(4*b*(a + b*x)*(-a*d + b*c)) - sqrt(c + d*x)/(2*b*(a + b*x)**2) + d**2*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(4*b**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1382():
    f = sqrt(c + d*x)/(a + b*x)**4
    F = d**2*sqrt(c + d*x)/(8*b*(a + b*x)*(-a*d + b*c)**2) - d*sqrt(c + d*x)/(12*b*(a + b*x)**2*(-a*d + b*c)) - sqrt(c + d*x)/(3*b*(a + b*x)**3) - d**3*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(8*b**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1383():
    f = sqrt(c + d*x)/(a + b*x)**5
    F = -5*d**3*sqrt(c + d*x)/(64*b*(a + b*x)*(-a*d + b*c)**3) + 5*d**2*sqrt(c + d*x)/(96*b*(a + b*x)**2*(-a*d + b*c)**2) - d*sqrt(c + d*x)/(24*b*(a + b*x)**3*(-a*d + b*c)) - sqrt(c + d*x)/(4*b*(a + b*x)**4) + 5*d**4*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(64*b**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1384():
    f = sqrt(c + d*x)/(a + b*x)**6
    F = 7*d**4*sqrt(c + d*x)/(128*b*(a + b*x)*(-a*d + b*c)**4) - 7*d**3*sqrt(c + d*x)/(192*b*(a + b*x)**2*(-a*d + b*c)**3) + 7*d**2*sqrt(c + d*x)/(240*b*(a + b*x)**3*(-a*d + b*c)**2) - d*sqrt(c + d*x)/(40*b*(a + b*x)**4*(-a*d + b*c)) - sqrt(c + d*x)/(5*b*(a + b*x)**5) - 7*d**5*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(128*b**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1385():
    f = (a + b*x)**5*(c + d*x)**(sympy.S(3)/2)
    F = 2*b**5*(c + d*x)**(sympy.S(15)/2)/(15*d**6) - 10*b**4*(c + d*x)**(sympy.S(13)/2)*(-a*d + b*c)/(13*d**6) + 20*b**3*(c + d*x)**(sympy.S(11)/2)*(-a*d + b*c)**2/(11*d**6) - 20*b**2*(c + d*x)**(sympy.S(9)/2)*(-a*d + b*c)**3/(9*d**6) + 10*b*(c + d*x)**(sympy.S(7)/2)*(-a*d + b*c)**4/(7*d**6) - 2*(c + d*x)**(sympy.S(5)/2)*(-a*d + b*c)**5/(5*d**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1386():
    f = (a + b*x)**4*(c + d*x)**(sympy.S(3)/2)
    F = 2*b**4*(c + d*x)**(sympy.S(13)/2)/(13*d**5) - 8*b**3*(c + d*x)**(sympy.S(11)/2)*(-a*d + b*c)/(11*d**5) + 4*b**2*(c + d*x)**(sympy.S(9)/2)*(-a*d + b*c)**2/(3*d**5) - 8*b*(c + d*x)**(sympy.S(7)/2)*(-a*d + b*c)**3/(7*d**5) + 2*(c + d*x)**(sympy.S(5)/2)*(-a*d + b*c)**4/(5*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1387():
    f = (a + b*x)**3*(c + d*x)**(sympy.S(3)/2)
    F = 2*b**3*(c + d*x)**(sympy.S(11)/2)/(11*d**4) - 2*b**2*(c + d*x)**(sympy.S(9)/2)*(-a*d + b*c)/(3*d**4) + 6*b*(c + d*x)**(sympy.S(7)/2)*(-a*d + b*c)**2/(7*d**4) - 2*(c + d*x)**(sympy.S(5)/2)*(-a*d + b*c)**3/(5*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1388():
    f = (a + b*x)**2*(c + d*x)**(sympy.S(3)/2)
    F = 2*b**2*(c + d*x)**(sympy.S(9)/2)/(9*d**3) - 4*b*(c + d*x)**(sympy.S(7)/2)*(-a*d + b*c)/(7*d**3) + 2*(c + d*x)**(sympy.S(5)/2)*(-a*d + b*c)**2/(5*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1389():
    f = (a + b*x)*(c + d*x)**(sympy.S(3)/2)
    F = 2*b*(c + d*x)**(sympy.S(7)/2)/(7*d**2) - (c + d*x)**(sympy.S(5)/2)*(-2*a*d + 2*b*c)/(5*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1390():
    f = (c + d*x)**(sympy.S(3)/2)
    F = 2*(c + d*x)**(sympy.S(5)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1391():
    f = (c + d*x)**(sympy.S(3)/2)/(a + b*x)
    F = 2*(c + d*x)**(sympy.S(3)/2)/(3*b) + sqrt(c + d*x)*(-2*a*d + 2*b*c)/b**2 - 2*(-a*d + b*c)**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1392():
    f = (c + d*x)**(sympy.S(3)/2)/(a + b*x)**2
    F = -(c + d*x)**(sympy.S(3)/2)/(b*(a + b*x)) + 3*d*sqrt(c + d*x)/b**2 - 3*d*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1393():
    f = (c + d*x)**(sympy.S(3)/2)/(a + b*x)**3
    F = -(c + d*x)**(sympy.S(3)/2)/(2*b*(a + b*x)**2) - 3*d*sqrt(c + d*x)/(4*b**2*(a + b*x)) - 3*d**2*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(4*b**(sympy.S(5)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1394():
    f = (c + d*x)**(sympy.S(3)/2)/(a + b*x)**4
    F = -(c + d*x)**(sympy.S(3)/2)/(3*b*(a + b*x)**3) - d**2*sqrt(c + d*x)/(8*b**2*(a + b*x)*(-a*d + b*c)) - d*sqrt(c + d*x)/(4*b**2*(a + b*x)**2) + d**3*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(8*b**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1395():
    f = (c + d*x)**(sympy.S(3)/2)/(a + b*x)**5
    F = -(c + d*x)**(sympy.S(3)/2)/(4*b*(a + b*x)**4) + 3*d**3*sqrt(c + d*x)/(64*b**2*(a + b*x)*(-a*d + b*c)**2) - d**2*sqrt(c + d*x)/(32*b**2*(a + b*x)**2*(-a*d + b*c)) - d*sqrt(c + d*x)/(8*b**2*(a + b*x)**3) - 3*d**4*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(64*b**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1396():
    f = (c + d*x)**(sympy.S(3)/2)/(a + b*x)**6
    F = -(c + d*x)**(sympy.S(3)/2)/(5*b*(a + b*x)**5) - 3*d**4*sqrt(c + d*x)/(128*b**2*(a + b*x)*(-a*d + b*c)**3) + d**3*sqrt(c + d*x)/(64*b**2*(a + b*x)**2*(-a*d + b*c)**2) - d**2*sqrt(c + d*x)/(80*b**2*(a + b*x)**3*(-a*d + b*c)) - 3*d*sqrt(c + d*x)/(40*b**2*(a + b*x)**4) + 3*d**5*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(128*b**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1397():
    f = (a + b*x)**5*(c + d*x)**(sympy.S(5)/2)
    F = 2*b**5*(c + d*x)**(sympy.S(17)/2)/(17*d**6) - 2*b**4*(c + d*x)**(sympy.S(15)/2)*(-a*d + b*c)/(3*d**6) + 20*b**3*(c + d*x)**(sympy.S(13)/2)*(-a*d + b*c)**2/(13*d**6) - 20*b**2*(c + d*x)**(sympy.S(11)/2)*(-a*d + b*c)**3/(11*d**6) + 10*b*(c + d*x)**(sympy.S(9)/2)*(-a*d + b*c)**4/(9*d**6) - 2*(c + d*x)**(sympy.S(7)/2)*(-a*d + b*c)**5/(7*d**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1398():
    f = (a + b*x)**4*(c + d*x)**(sympy.S(5)/2)
    F = 2*b**4*(c + d*x)**(sympy.S(15)/2)/(15*d**5) - 8*b**3*(c + d*x)**(sympy.S(13)/2)*(-a*d + b*c)/(13*d**5) + 12*b**2*(c + d*x)**(sympy.S(11)/2)*(-a*d + b*c)**2/(11*d**5) - 8*b*(c + d*x)**(sympy.S(9)/2)*(-a*d + b*c)**3/(9*d**5) + 2*(c + d*x)**(sympy.S(7)/2)*(-a*d + b*c)**4/(7*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1399():
    f = (a + b*x)**3*(c + d*x)**(sympy.S(5)/2)
    F = 2*b**3*(c + d*x)**(sympy.S(13)/2)/(13*d**4) - 6*b**2*(c + d*x)**(sympy.S(11)/2)*(-a*d + b*c)/(11*d**4) + 2*b*(c + d*x)**(sympy.S(9)/2)*(-a*d + b*c)**2/(3*d**4) - 2*(c + d*x)**(sympy.S(7)/2)*(-a*d + b*c)**3/(7*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1400():
    f = (a + b*x)**2*(c + d*x)**(sympy.S(5)/2)
    F = 2*b**2*(c + d*x)**(sympy.S(11)/2)/(11*d**3) - 4*b*(c + d*x)**(sympy.S(9)/2)*(-a*d + b*c)/(9*d**3) + 2*(c + d*x)**(sympy.S(7)/2)*(-a*d + b*c)**2/(7*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1401():
    f = (a + b*x)*(c + d*x)**(sympy.S(5)/2)
    F = 2*b*(c + d*x)**(sympy.S(9)/2)/(9*d**2) - (c + d*x)**(sympy.S(7)/2)*(-2*a*d + 2*b*c)/(7*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1402():
    f = (c + d*x)**(sympy.S(5)/2)
    F = 2*(c + d*x)**(sympy.S(7)/2)/(7*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1403():
    f = (c + d*x)**(sympy.S(5)/2)/(a + b*x)
    F = 2*(c + d*x)**(sympy.S(5)/2)/(5*b) + (c + d*x)**(sympy.S(3)/2)*(-2*a*d + 2*b*c)/(3*b**2) + 2*sqrt(c + d*x)*(-a*d + b*c)**2/b**3 - 2*(-a*d + b*c)**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1404():
    f = (c + d*x)**(sympy.S(5)/2)/(a + b*x)**2
    F = -(c + d*x)**(sympy.S(5)/2)/(b*(a + b*x)) + 5*d*(c + d*x)**(sympy.S(3)/2)/(3*b**2) + 5*d*sqrt(c + d*x)*(-a*d + b*c)/b**3 - 5*d*(-a*d + b*c)**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1405():
    f = (c + d*x)**(sympy.S(5)/2)/(a + b*x)**3
    F = -(c + d*x)**(sympy.S(5)/2)/(2*b*(a + b*x)**2) - 5*d*(c + d*x)**(sympy.S(3)/2)/(4*b**2*(a + b*x)) + 15*d**2*sqrt(c + d*x)/(4*b**3) - 15*d**2*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(4*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1406():
    f = (c + d*x)**(sympy.S(5)/2)/(a + b*x)**4
    F = -(c + d*x)**(sympy.S(5)/2)/(3*b*(a + b*x)**3) - 5*d*(c + d*x)**(sympy.S(3)/2)/(12*b**2*(a + b*x)**2) - 5*d**2*sqrt(c + d*x)/(8*b**3*(a + b*x)) - 5*d**3*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(8*b**(sympy.S(7)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1407():
    f = (c + d*x)**(sympy.S(5)/2)/(a + b*x)**5
    F = -(c + d*x)**(sympy.S(5)/2)/(4*b*(a + b*x)**4) - 5*d*(c + d*x)**(sympy.S(3)/2)/(24*b**2*(a + b*x)**3) - 5*d**3*sqrt(c + d*x)/(64*b**3*(a + b*x)*(-a*d + b*c)) - 5*d**2*sqrt(c + d*x)/(32*b**3*(a + b*x)**2) + 5*d**4*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(64*b**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1408():
    f = (c + d*x)**(sympy.S(5)/2)/(a + b*x)**6
    F = -(c + d*x)**(sympy.S(5)/2)/(5*b*(a + b*x)**5) - d*(c + d*x)**(sympy.S(3)/2)/(8*b**2*(a + b*x)**4) + 3*d**4*sqrt(c + d*x)/(128*b**3*(a + b*x)*(-a*d + b*c)**2) - d**3*sqrt(c + d*x)/(64*b**3*(a + b*x)**2*(-a*d + b*c)) - d**2*sqrt(c + d*x)/(16*b**3*(a + b*x)**3) - 3*d**5*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(128*b**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1409():
    f = sqrt(x - 1)/(x + 1)**2
    F = -sqrt(x - 1)/(x + 1) + sqrt(2)*atan(sqrt(2)*sqrt(x - 1)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1410():
    f = sqrt(x - 1)/(x + 1)**3
    F = sqrt(x - 1)/(8*x + 8) - sqrt(x - 1)/(2*(x + 1)**2) + sqrt(2)*atan(sqrt(2)*sqrt(x - 1)/2)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1411():
    f = (a + b*x)**5/sqrt(c + d*x)
    F = 2*b**5*(c + d*x)**(sympy.S(11)/2)/(11*d**6) - 10*b**4*(c + d*x)**(sympy.S(9)/2)*(-a*d + b*c)/(9*d**6) + 20*b**3*(c + d*x)**(sympy.S(7)/2)*(-a*d + b*c)**2/(7*d**6) - 4*b**2*(c + d*x)**(sympy.S(5)/2)*(-a*d + b*c)**3/d**6 + 10*b*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**4/(3*d**6) - 2*sqrt(c + d*x)*(-a*d + b*c)**5/d**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1412():
    f = (a + b*x)**4/sqrt(c + d*x)
    F = 2*b**4*(c + d*x)**(sympy.S(9)/2)/(9*d**5) - 8*b**3*(c + d*x)**(sympy.S(7)/2)*(-a*d + b*c)/(7*d**5) + 12*b**2*(c + d*x)**(sympy.S(5)/2)*(-a*d + b*c)**2/(5*d**5) - 8*b*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**3/(3*d**5) + 2*sqrt(c + d*x)*(-a*d + b*c)**4/d**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1413():
    f = (a + b*x)**3/sqrt(c + d*x)
    F = 2*b**3*(c + d*x)**(sympy.S(7)/2)/(7*d**4) - 6*b**2*(c + d*x)**(sympy.S(5)/2)*(-a*d + b*c)/(5*d**4) + 2*b*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**2/d**4 - 2*sqrt(c + d*x)*(-a*d + b*c)**3/d**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1414():
    f = (a + b*x)**2/sqrt(c + d*x)
    F = 2*b**2*(c + d*x)**(sympy.S(5)/2)/(5*d**3) - 4*b*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)/(3*d**3) + 2*sqrt(c + d*x)*(-a*d + b*c)**2/d**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1415():
    f = (a + b*x)/sqrt(c + d*x)
    F = 2*b*(c + d*x)**(sympy.S(3)/2)/(3*d**2) - sqrt(c + d*x)*(-2*a*d + 2*b*c)/d**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1416():
    f = 1/sqrt(c + d*x)
    F = 2*sqrt(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1417():
    f = 1/((a + b*x)*sqrt(c + d*x))
    F = -2*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(sqrt(b)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1418():
    f = 1/((a + b*x)**2*sqrt(c + d*x))
    F = -sqrt(c + d*x)/((a + b*x)*(-a*d + b*c)) + d*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(sqrt(b)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1419():
    f = 1/((a + b*x)**3*sqrt(c + d*x))
    F = 3*d*sqrt(c + d*x)/(4*(a + b*x)*(-a*d + b*c)**2) - sqrt(c + d*x)/((a + b*x)**2*(-2*a*d + 2*b*c)) - 3*d**2*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(4*sqrt(b)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1420():
    f = 1/((a + b*x)**4*sqrt(c + d*x))
    F = -5*d**2*sqrt(c + d*x)/(8*(a + b*x)*(-a*d + b*c)**3) + 5*d*sqrt(c + d*x)/(12*(a + b*x)**2*(-a*d + b*c)**2) - sqrt(c + d*x)/((a + b*x)**3*(-3*a*d + 3*b*c)) + 5*d**3*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(8*sqrt(b)*(-a*d + b*c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1421():
    f = 1/((a + b*x)**5*sqrt(c + d*x))
    F = 35*d**3*sqrt(c + d*x)/(64*(a + b*x)*(-a*d + b*c)**4) - 35*d**2*sqrt(c + d*x)/(96*(a + b*x)**2*(-a*d + b*c)**3) + 7*d*sqrt(c + d*x)/(24*(a + b*x)**3*(-a*d + b*c)**2) - sqrt(c + d*x)/((a + b*x)**4*(-4*a*d + 4*b*c)) - 35*d**4*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(64*sqrt(b)*(-a*d + b*c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1422():
    f = (a + b*x)**5/(c + d*x)**(sympy.S(3)/2)
    F = 2*b**5*(c + d*x)**(sympy.S(9)/2)/(9*d**6) - 10*b**4*(c + d*x)**(sympy.S(7)/2)*(-a*d + b*c)/(7*d**6) + 4*b**3*(c + d*x)**(sympy.S(5)/2)*(-a*d + b*c)**2/d**6 - 20*b**2*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**3/(3*d**6) + 10*b*sqrt(c + d*x)*(-a*d + b*c)**4/d**6 + 2*(-a*d + b*c)**5/(d**6*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1423():
    f = (a + b*x)**4/(c + d*x)**(sympy.S(3)/2)
    F = 2*b**4*(c + d*x)**(sympy.S(7)/2)/(7*d**5) - 8*b**3*(c + d*x)**(sympy.S(5)/2)*(-a*d + b*c)/(5*d**5) + 4*b**2*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**2/d**5 - 8*b*sqrt(c + d*x)*(-a*d + b*c)**3/d**5 - 2*(-a*d + b*c)**4/(d**5*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1424():
    f = (a + b*x)**3/(c + d*x)**(sympy.S(3)/2)
    F = 2*b**3*(c + d*x)**(sympy.S(5)/2)/(5*d**4) - 2*b**2*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)/d**4 + 6*b*sqrt(c + d*x)*(-a*d + b*c)**2/d**4 + 2*(-a*d + b*c)**3/(d**4*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1425():
    f = (a + b*x)**2/(c + d*x)**(sympy.S(3)/2)
    F = 2*b**2*(c + d*x)**(sympy.S(3)/2)/(3*d**3) - 4*b*sqrt(c + d*x)*(-a*d + b*c)/d**3 - 2*(-a*d + b*c)**2/(d**3*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1426():
    f = (a + b*x)/(c + d*x)**(sympy.S(3)/2)
    F = 2*b*sqrt(c + d*x)/d**2 + (-2*a*d + 2*b*c)/(d**2*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1427():
    f = (c + d*x)**(sympy.S(-3)/2)
    F = -2/(d*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1428():
    f = 1/((a + b*x)*(c + d*x)**(sympy.S(3)/2))
    F = -2*sqrt(b)*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(-a*d + b*c)**(sympy.S(3)/2) + 2/(sqrt(c + d*x)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1429():
    f = 1/((a + b*x)**2*(c + d*x)**(sympy.S(3)/2))
    F = 3*sqrt(b)*d*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(-a*d + b*c)**(sympy.S(5)/2) - 3*d/(sqrt(c + d*x)*(-a*d + b*c)**2) - 1/((a + b*x)*sqrt(c + d*x)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1430():
    f = 1/((a + b*x)**3*(c + d*x)**(sympy.S(3)/2))
    F = -15*sqrt(b)*d**2*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(4*(-a*d + b*c)**(sympy.S(7)/2)) + 15*d**2/(4*sqrt(c + d*x)*(-a*d + b*c)**3) + 5*d/(4*(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**2) - 1/((a + b*x)**2*sqrt(c + d*x)*(-2*a*d + 2*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1431():
    f = 1/((a + b*x)**4*(c + d*x)**(sympy.S(3)/2))
    F = 35*sqrt(b)*d**3*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(8*(-a*d + b*c)**(sympy.S(9)/2)) - 35*d**3/(8*sqrt(c + d*x)*(-a*d + b*c)**4) - 35*d**2/(24*(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**3) + 7*d/(12*(a + b*x)**2*sqrt(c + d*x)*(-a*d + b*c)**2) - 1/((a + b*x)**3*sqrt(c + d*x)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1432():
    f = (a + b*x)**5/(c + d*x)**(sympy.S(5)/2)
    F = 2*b**5*(c + d*x)**(sympy.S(7)/2)/(7*d**6) - 2*b**4*(c + d*x)**(sympy.S(5)/2)*(-a*d + b*c)/d**6 + 20*b**3*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**2/(3*d**6) - 20*b**2*sqrt(c + d*x)*(-a*d + b*c)**3/d**6 - 10*b*(-a*d + b*c)**4/(d**6*sqrt(c + d*x)) + 2*(-a*d + b*c)**5/(3*d**6*(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1433():
    f = (a + b*x)**4/(c + d*x)**(sympy.S(5)/2)
    F = 2*b**4*(c + d*x)**(sympy.S(5)/2)/(5*d**5) - 8*b**3*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)/(3*d**5) + 12*b**2*sqrt(c + d*x)*(-a*d + b*c)**2/d**5 + 8*b*(-a*d + b*c)**3/(d**5*sqrt(c + d*x)) - 2*(-a*d + b*c)**4/(3*d**5*(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1434():
    f = (a + b*x)**3/(c + d*x)**(sympy.S(5)/2)
    F = 2*b**3*(c + d*x)**(sympy.S(3)/2)/(3*d**4) - 6*b**2*sqrt(c + d*x)*(-a*d + b*c)/d**4 - 6*b*(-a*d + b*c)**2/(d**4*sqrt(c + d*x)) + 2*(-a*d + b*c)**3/(3*d**4*(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1435():
    f = (a + b*x)**2/(c + d*x)**(sympy.S(5)/2)
    F = 2*b**2*sqrt(c + d*x)/d**3 + 4*b*(-a*d + b*c)/(d**3*sqrt(c + d*x)) - 2*(-a*d + b*c)**2/(3*d**3*(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1436():
    f = (a + b*x)/(c + d*x)**(sympy.S(5)/2)
    F = -2*b/(d**2*sqrt(c + d*x)) + (-2*a*d + 2*b*c)/(3*d**2*(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1437():
    f = (c + d*x)**(sympy.S(-5)/2)
    F = -2/(3*d*(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1438():
    f = 1/((a + b*x)*(c + d*x)**(sympy.S(5)/2))
    F = -2*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(-a*d + b*c)**(sympy.S(5)/2) + 2*b/(sqrt(c + d*x)*(-a*d + b*c)**2) + 2/((c + d*x)**(sympy.S(3)/2)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1439():
    f = 1/((a + b*x)**2*(c + d*x)**(sympy.S(5)/2))
    F = 5*b**(sympy.S(3)/2)*d*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(-a*d + b*c)**(sympy.S(7)/2) - 5*b*d/(sqrt(c + d*x)*(-a*d + b*c)**3) - 5*d/(3*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**2) - 1/((a + b*x)*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1440():
    f = 1/((a + b*x)**3*(c + d*x)**(sympy.S(5)/2))
    F = -35*b**(sympy.S(3)/2)*d**2*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(4*(-a*d + b*c)**(sympy.S(9)/2)) + 35*b*d**2/(4*sqrt(c + d*x)*(-a*d + b*c)**4) + 35*d**2/(12*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**3) + 7*d/(4*(a + b*x)*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**2) - 1/((a + b*x)**2*(c + d*x)**(sympy.S(3)/2)*(-2*a*d + 2*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1441():
    f = 1/((a + b*x)**4*(c + d*x)**(sympy.S(5)/2))
    F = 105*b**(sympy.S(3)/2)*d**3*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(8*(-a*d + b*c)**(sympy.S(11)/2)) - 105*b*d**3/(8*sqrt(c + d*x)*(-a*d + b*c)**5) - 35*d**3/(8*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**4) - 21*d**2/(8*(a + b*x)*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**3) + 3*d/(4*(a + b*x)**2*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**2) - 1/((a + b*x)**3*(c + d*x)**(sympy.S(3)/2)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1442():
    f = (a + b*x)**5*(a*c + b*c*x)**(sympy.S(3)/2)
    F = 2*(a*c + b*c*x)**(sympy.S(15)/2)/(15*b*c**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1443():
    f = (a + b*x)**5*sqrt(a*c + b*c*x)
    F = 2*(a*c + b*c*x)**(sympy.S(13)/2)/(13*b*c**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1444():
    f = (a + b*x)**5/sqrt(a*c + b*c*x)
    F = 2*(a*c + b*c*x)**(sympy.S(11)/2)/(11*b*c**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1445():
    f = (a + b*x)**5/(a*c + b*c*x)**(sympy.S(3)/2)
    F = 2*(a*c + b*c*x)**(sympy.S(9)/2)/(9*b*c**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1446():
    f = (a + b*x)**5/(a*c + b*c*x)**(sympy.S(5)/2)
    F = 2*(a*c + b*c*x)**(sympy.S(7)/2)/(7*b*c**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1447():
    f = (a + b*x)**5/(a*c + b*c*x)**(sympy.S(7)/2)
    F = 2*(a*c + b*c*x)**(sympy.S(5)/2)/(5*b*c**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1448():
    f = (a + b*x)**5/(a*c + b*c*x)**(sympy.S(9)/2)
    F = 2*(a*c + b*c*x)**(sympy.S(3)/2)/(3*b*c**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1449():
    f = (a + b*x)**5/(a*c + b*c*x)**(sympy.S(11)/2)
    F = 2*sqrt(a*c + b*c*x)/(b*c**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1450():
    f = (a + b*x)**5/(a*c + b*c*x)**(sympy.S(13)/2)
    F = -2/(b*c**6*sqrt(a*c + b*c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1451():
    f = 1/((x - 2)*sqrt(x + 2))
    F = -atanh(sqrt(x + 2)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1452():
    f = 1/((3*x + 2)*sqrt(5*x + 1))
    F = 2*sqrt(21)*atan(sqrt(21)*sqrt(5*x + 1)/7)/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1453():
    f = (1 - x)**(sympy.S(1)/3)/(x + 1)
    F = 3*(1 - x)**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*log(x + 1)/2 + 3*2**(sympy.S(1)/3)*log(-(1 - x)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/2 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x)**(sympy.S(1)/3) + 1)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1454():
    f = (3 - 2*x)**(sympy.S(1)/3)*(x + 7)
    F = 3*(3 - 2*x)**(sympy.S(7)/3)/28 - 51*(3 - 2*x)**(sympy.S(4)/3)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1455():
    f = (1 - x)**(sympy.S(1)/3)*(x + 1)**2
    F = -3*(1 - x)**(sympy.S(10)/3)/10 + 12*(1 - x)**(sympy.S(7)/3)/7 - 3*(1 - x)**(sympy.S(4)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1456():
    f = 1/((a + b*x)*(c + d*x)**(sympy.S(1)/3))
    F = -log(a + b*x)/(2*b**(sympy.S(2)/3)*(-a*d + b*c)**(sympy.S(1)/3)) + 3*log(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(2*b**(sympy.S(2)/3)*(-a*d + b*c)**(sympy.S(1)/3)) + sqrt(3)*atan(sqrt(3)*(2*b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(-a*d + b*c)**(sympy.S(1)/3) + 1)/3)/(b**(sympy.S(2)/3)*(-a*d + b*c)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1457():
    f = 1/((a + b*x)*(c + d*x)**(sympy.S(2)/3))
    F = -log(a + b*x)/(2*b**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3)) + 3*log(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(2*b**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3)) - sqrt(3)*atan(sqrt(3)*(2*b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(-a*d + b*c)**(sympy.S(1)/3) + 1)/3)/(b**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1458():
    f = (a + b*x)**(sympy.S(7)/2)*sqrt(c + d*x)
    F = (a + b*x)**(sympy.S(9)/2)*sqrt(c + d*x)/(5*b) + (a + b*x)**(sympy.S(7)/2)*sqrt(c + d*x)*(-a*d + b*c)/(40*b*d) - 7*(a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)*(-a*d + b*c)**2/(240*b*d**2) + 7*(a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-a*d + b*c)**3/(192*b*d**3) - 7*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**4/(128*b*d**4) + 7*(-a*d + b*c)**5*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(128*b**(sympy.S(3)/2)*d**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1459():
    f = (a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)
    F = (a + b*x)**(sympy.S(7)/2)*sqrt(c + d*x)/(4*b) + (a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)*(-a*d + b*c)/(24*b*d) - 5*(a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-a*d + b*c)**2/(96*b*d**2) + 5*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**3/(64*b*d**3) - 5*(-a*d + b*c)**4*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(64*b**(sympy.S(3)/2)*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1460():
    f = (a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)
    F = (a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)/(3*b) + (a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-a*d + b*c)/(12*b*d) - sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**2/(8*b*d**2) + (-a*d + b*c)**3*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(8*b**(sympy.S(3)/2)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1461():
    f = sqrt(a + b*x)*sqrt(c + d*x)
    F = (a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)/(2*b) + sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)/(4*b*d) - (-a*d + b*c)**2*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(4*b**(sympy.S(3)/2)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1462():
    f = sqrt(c + d*x)/sqrt(a + b*x)
    F = sqrt(a + b*x)*sqrt(c + d*x)/b + (-a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(b**(sympy.S(3)/2)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1463():
    f = sqrt(c + d*x)/(a + b*x)**(sympy.S(3)/2)
    F = -2*sqrt(c + d*x)/(b*sqrt(a + b*x)) + 2*sqrt(d)*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1464():
    f = sqrt(c + d*x)/(a + b*x)**(sympy.S(5)/2)
    F = -2*(c + d*x)**(sympy.S(3)/2)/((a + b*x)**(sympy.S(3)/2)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1465():
    f = sqrt(c + d*x)/(a + b*x)**(sympy.S(7)/2)
    F = 4*d*(c + d*x)**(sympy.S(3)/2)/(15*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(3)/2)/((a + b*x)**(sympy.S(5)/2)*(-5*a*d + 5*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1466():
    f = sqrt(c + d*x)/(a + b*x)**(sympy.S(9)/2)
    F = -16*d**2*(c + d*x)**(sympy.S(3)/2)/(105*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)**3) + 8*d*(c + d*x)**(sympy.S(3)/2)/(35*(a + b*x)**(sympy.S(5)/2)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(3)/2)/((a + b*x)**(sympy.S(7)/2)*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1467():
    f = sqrt(c + d*x)/(a + b*x)**(sympy.S(11)/2)
    F = 32*d**3*(c + d*x)**(sympy.S(3)/2)/(315*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)**4) - 16*d**2*(c + d*x)**(sympy.S(3)/2)/(105*(a + b*x)**(sympy.S(5)/2)*(-a*d + b*c)**3) + 4*d*(c + d*x)**(sympy.S(3)/2)/(21*(a + b*x)**(sympy.S(7)/2)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(3)/2)/((a + b*x)**(sympy.S(9)/2)*(-9*a*d + 9*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1468():
    f = sqrt(c + d*x)/(a + b*x)**(sympy.S(13)/2)
    F = -256*d**4*(c + d*x)**(sympy.S(3)/2)/(3465*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)**5) + 128*d**3*(c + d*x)**(sympy.S(3)/2)/(1155*(a + b*x)**(sympy.S(5)/2)*(-a*d + b*c)**4) - 32*d**2*(c + d*x)**(sympy.S(3)/2)/(231*(a + b*x)**(sympy.S(7)/2)*(-a*d + b*c)**3) + 16*d*(c + d*x)**(sympy.S(3)/2)/(99*(a + b*x)**(sympy.S(9)/2)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(3)/2)/((a + b*x)**(sympy.S(11)/2)*(-11*a*d + 11*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1469():
    f = (a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(3)/2)
    F = (a + b*x)**(sympy.S(7)/2)*(c + d*x)**(sympy.S(3)/2)/(5*b) + (a + b*x)**(sympy.S(7)/2)*sqrt(c + d*x)*(-3*a*d + 3*b*c)/(40*b**2) + (a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)*(-a*d + b*c)**2/(80*b**2*d) - (a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-a*d + b*c)**3/(64*b**2*d**2) + 3*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**4/(128*b**2*d**3) - 3*(-a*d + b*c)**5*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(128*b**(sympy.S(5)/2)*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1470():
    f = (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)
    F = (a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(3)/2)/(4*b) + (a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)*(-a*d + b*c)/(8*b**2) + (a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-a*d + b*c)**2/(32*b**2*d) - 3*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**3/(64*b**2*d**2) + 3*(-a*d + b*c)**4*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(64*b**(sympy.S(5)/2)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1471():
    f = sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/2)
    F = (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)/(3*b) + (a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-a*d + b*c)/(4*b**2) + sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**2/(8*b**2*d) - (-a*d + b*c)**3*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(8*b**(sympy.S(5)/2)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1472():
    f = (c + d*x)**(sympy.S(3)/2)/sqrt(a + b*x)
    F = sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/2)/(2*b) + sqrt(a + b*x)*sqrt(c + d*x)*(-3*a*d + 3*b*c)/(4*b**2) + 3*(-a*d + b*c)**2*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(4*b**(sympy.S(5)/2)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1473():
    f = (c + d*x)**(sympy.S(3)/2)/(a + b*x)**(sympy.S(3)/2)
    F = -2*(c + d*x)**(sympy.S(3)/2)/(b*sqrt(a + b*x)) + 3*d*sqrt(a + b*x)*sqrt(c + d*x)/b**2 + 3*sqrt(d)*(-a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1474():
    f = (c + d*x)**(sympy.S(3)/2)/(a + b*x)**(sympy.S(5)/2)
    F = -2*(c + d*x)**(sympy.S(3)/2)/(3*b*(a + b*x)**(sympy.S(3)/2)) - 2*d*sqrt(c + d*x)/(b**2*sqrt(a + b*x)) + 2*d**(sympy.S(3)/2)*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1475():
    f = (c + d*x)**(sympy.S(3)/2)/(a + b*x)**(sympy.S(7)/2)
    F = -2*(c + d*x)**(sympy.S(5)/2)/((a + b*x)**(sympy.S(5)/2)*(-5*a*d + 5*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1476():
    f = (c + d*x)**(sympy.S(3)/2)/(a + b*x)**(sympy.S(9)/2)
    F = 4*d*(c + d*x)**(sympy.S(5)/2)/(35*(a + b*x)**(sympy.S(5)/2)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(5)/2)/((a + b*x)**(sympy.S(7)/2)*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1477():
    f = (c + d*x)**(sympy.S(3)/2)/(a + b*x)**(sympy.S(11)/2)
    F = -16*d**2*(c + d*x)**(sympy.S(5)/2)/(315*(a + b*x)**(sympy.S(5)/2)*(-a*d + b*c)**3) + 8*d*(c + d*x)**(sympy.S(5)/2)/(63*(a + b*x)**(sympy.S(7)/2)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(5)/2)/((a + b*x)**(sympy.S(9)/2)*(-9*a*d + 9*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1478():
    f = (c + d*x)**(sympy.S(3)/2)/(a + b*x)**(sympy.S(13)/2)
    F = 32*d**3*(c + d*x)**(sympy.S(5)/2)/(1155*(a + b*x)**(sympy.S(5)/2)*(-a*d + b*c)**4) - 16*d**2*(c + d*x)**(sympy.S(5)/2)/(231*(a + b*x)**(sympy.S(7)/2)*(-a*d + b*c)**3) + 4*d*(c + d*x)**(sympy.S(5)/2)/(33*(a + b*x)**(sympy.S(9)/2)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(5)/2)/((a + b*x)**(sympy.S(11)/2)*(-11*a*d + 11*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1479():
    f = (a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(5)/2)
    F = (a + b*x)**(sympy.S(7)/2)*(c + d*x)**(sympy.S(5)/2)/(6*b) + (a + b*x)**(sympy.S(7)/2)*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)/(12*b**2) + (a + b*x)**(sympy.S(7)/2)*sqrt(c + d*x)*(-a*d + b*c)**2/(32*b**3) + (a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)*(-a*d + b*c)**3/(192*b**3*d) - 5*(a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-a*d + b*c)**4/(768*b**3*d**2) + 5*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**5/(512*b**3*d**3) - 5*(-a*d + b*c)**6*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(512*b**(sympy.S(7)/2)*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1480():
    f = (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(5)/2)
    F = (a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(5)/2)/(5*b) + (a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)/(8*b**2) + (a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)*(-a*d + b*c)**2/(16*b**3) + (a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-a*d + b*c)**3/(64*b**3*d) - 3*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**4/(128*b**3*d**2) + 3*(-a*d + b*c)**5*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(128*b**(sympy.S(7)/2)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1481():
    f = sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/2)
    F = (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(5)/2)/(4*b) + (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)*(-5*a*d + 5*b*c)/(24*b**2) + 5*(a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-a*d + b*c)**2/(32*b**3) + 5*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**3/(64*b**3*d) - 5*(-a*d + b*c)**4*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(64*b**(sympy.S(7)/2)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1482():
    f = (c + d*x)**(sympy.S(5)/2)/sqrt(a + b*x)
    F = sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/2)/(3*b) + sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/2)*(-5*a*d + 5*b*c)/(12*b**2) + 5*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**2/(8*b**3) + 5*(-a*d + b*c)**3*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(8*b**(sympy.S(7)/2)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1483():
    f = (c + d*x)**(sympy.S(5)/2)/(a + b*x)**(sympy.S(3)/2)
    F = -2*(c + d*x)**(sympy.S(5)/2)/(b*sqrt(a + b*x)) + 5*d*sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/2)/(2*b**2) + 15*d*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)/(4*b**3) + 15*sqrt(d)*(-a*d + b*c)**2*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(4*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1484():
    f = (c + d*x)**(sympy.S(5)/2)/(a + b*x)**(sympy.S(5)/2)
    F = -2*(c + d*x)**(sympy.S(5)/2)/(3*b*(a + b*x)**(sympy.S(3)/2)) - 10*d*(c + d*x)**(sympy.S(3)/2)/(3*b**2*sqrt(a + b*x)) + 5*d**2*sqrt(a + b*x)*sqrt(c + d*x)/b**3 + 5*d**(sympy.S(3)/2)*(-a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1485():
    f = (c + d*x)**(sympy.S(5)/2)/(a + b*x)**(sympy.S(7)/2)
    F = -2*(c + d*x)**(sympy.S(5)/2)/(5*b*(a + b*x)**(sympy.S(5)/2)) - 2*d*(c + d*x)**(sympy.S(3)/2)/(3*b**2*(a + b*x)**(sympy.S(3)/2)) - 2*d**2*sqrt(c + d*x)/(b**3*sqrt(a + b*x)) + 2*d**(sympy.S(5)/2)*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1486():
    f = (c + d*x)**(sympy.S(5)/2)/(a + b*x)**(sympy.S(9)/2)
    F = -2*(c + d*x)**(sympy.S(7)/2)/((a + b*x)**(sympy.S(7)/2)*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1487():
    f = (c + d*x)**(sympy.S(5)/2)/(a + b*x)**(sympy.S(11)/2)
    F = 4*d*(c + d*x)**(sympy.S(7)/2)/(63*(a + b*x)**(sympy.S(7)/2)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(7)/2)/((a + b*x)**(sympy.S(9)/2)*(-9*a*d + 9*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1488():
    f = (c + d*x)**(sympy.S(5)/2)/(a + b*x)**(sympy.S(13)/2)
    F = -16*d**2*(c + d*x)**(sympy.S(7)/2)/(693*(a + b*x)**(sympy.S(7)/2)*(-a*d + b*c)**3) + 8*d*(c + d*x)**(sympy.S(7)/2)/(99*(a + b*x)**(sympy.S(9)/2)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(7)/2)/((a + b*x)**(sympy.S(11)/2)*(-11*a*d + 11*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1489():
    f = (c + d*x)**(sympy.S(5)/2)/(a + b*x)**(sympy.S(15)/2)
    F = 32*d**3*(c + d*x)**(sympy.S(7)/2)/(3003*(a + b*x)**(sympy.S(7)/2)*(-a*d + b*c)**4) - 16*d**2*(c + d*x)**(sympy.S(7)/2)/(429*(a + b*x)**(sympy.S(9)/2)*(-a*d + b*c)**3) + 12*d*(c + d*x)**(sympy.S(7)/2)/(143*(a + b*x)**(sympy.S(11)/2)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(7)/2)/((a + b*x)**(sympy.S(13)/2)*(-13*a*d + 13*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1490():
    f = (a + b*x)**(sympy.S(7)/2)/sqrt(c + d*x)
    F = (a + b*x)**(sympy.S(7)/2)*sqrt(c + d*x)/(4*d) - (a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)*(-7*a*d + 7*b*c)/(24*d**2) + 35*(a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-a*d + b*c)**2/(96*d**3) - 35*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**3/(64*d**4) + 35*(-a*d + b*c)**4*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(64*sqrt(b)*d**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1491():
    f = (a + b*x)**(sympy.S(5)/2)/sqrt(c + d*x)
    F = (a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)/(3*d) - (a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-5*a*d + 5*b*c)/(12*d**2) + 5*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**2/(8*d**3) - 5*(-a*d + b*c)**3*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(8*sqrt(b)*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1492():
    f = (a + b*x)**(sympy.S(3)/2)/sqrt(c + d*x)
    F = (a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)/(2*d) - sqrt(a + b*x)*sqrt(c + d*x)*(-3*a*d + 3*b*c)/(4*d**2) + 3*(-a*d + b*c)**2*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(4*sqrt(b)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1493():
    f = sqrt(a + b*x)/sqrt(c + d*x)
    F = sqrt(a + b*x)*sqrt(c + d*x)/d - (-a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(sqrt(b)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1494():
    f = 1/(sqrt(a + b*x)*sqrt(c + d*x))
    F = 2*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(sqrt(b)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1495():
    f = 1/((a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x))
    F = -2*sqrt(c + d*x)/(sqrt(a + b*x)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1496():
    f = 1/((a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x))
    F = 4*d*sqrt(c + d*x)/(3*sqrt(a + b*x)*(-a*d + b*c)**2) - 2*sqrt(c + d*x)/((a + b*x)**(sympy.S(3)/2)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1497():
    f = 1/((a + b*x)**(sympy.S(7)/2)*sqrt(c + d*x))
    F = -16*d**2*sqrt(c + d*x)/(15*sqrt(a + b*x)*(-a*d + b*c)**3) + 8*d*sqrt(c + d*x)/(15*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)**2) - 2*sqrt(c + d*x)/((a + b*x)**(sympy.S(5)/2)*(-5*a*d + 5*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1498():
    f = 1/((a + b*x)**(sympy.S(9)/2)*sqrt(c + d*x))
    F = 32*d**3*sqrt(c + d*x)/(35*sqrt(a + b*x)*(-a*d + b*c)**4) - 16*d**2*sqrt(c + d*x)/(35*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)**3) + 12*d*sqrt(c + d*x)/(35*(a + b*x)**(sympy.S(5)/2)*(-a*d + b*c)**2) - 2*sqrt(c + d*x)/((a + b*x)**(sympy.S(7)/2)*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1499():
    f = 1/((a + b*x)**(sympy.S(11)/2)*sqrt(c + d*x))
    F = -256*d**4*sqrt(c + d*x)/(315*sqrt(a + b*x)*(-a*d + b*c)**5) + 128*d**3*sqrt(c + d*x)/(315*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)**4) - 32*d**2*sqrt(c + d*x)/(105*(a + b*x)**(sympy.S(5)/2)*(-a*d + b*c)**3) + 16*d*sqrt(c + d*x)/(63*(a + b*x)**(sympy.S(7)/2)*(-a*d + b*c)**2) - 2*sqrt(c + d*x)/((a + b*x)**(sympy.S(9)/2)*(-9*a*d + 9*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1500():
    f = (a + b*x)**(sympy.S(7)/2)/(c + d*x)**(sympy.S(3)/2)
    F = -35*sqrt(b)*(-a*d + b*c)**3*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(8*d**(sympy.S(9)/2)) + 7*b*(a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)/(3*d**2) - 35*b*(a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-a*d + b*c)/(12*d**3) + 35*b*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**2/(8*d**4) - 2*(a + b*x)**(sympy.S(7)/2)/(d*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1501():
    f = (a + b*x)**(sympy.S(5)/2)/(c + d*x)**(sympy.S(3)/2)
    F = 15*sqrt(b)*(-a*d + b*c)**2*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(4*d**(sympy.S(7)/2)) + 5*b*(a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)/(2*d**2) - 15*b*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)/(4*d**3) - 2*(a + b*x)**(sympy.S(5)/2)/(d*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1502():
    f = (a + b*x)**(sympy.S(3)/2)/(c + d*x)**(sympy.S(3)/2)
    F = -3*sqrt(b)*(-a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/d**(sympy.S(5)/2) + 3*b*sqrt(a + b*x)*sqrt(c + d*x)/d**2 - 2*(a + b*x)**(sympy.S(3)/2)/(d*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1503():
    f = sqrt(a + b*x)/(c + d*x)**(sympy.S(3)/2)
    F = 2*sqrt(b)*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/d**(sympy.S(3)/2) - 2*sqrt(a + b*x)/(d*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1504():
    f = 1/(sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/2))
    F = 2*sqrt(a + b*x)/(sqrt(c + d*x)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1505():
    f = 1/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2))
    F = -4*d*sqrt(a + b*x)/(sqrt(c + d*x)*(-a*d + b*c)**2) - 2/(sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1506():
    f = 1/((a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(3)/2))
    F = 16*d**2*sqrt(a + b*x)/(3*sqrt(c + d*x)*(-a*d + b*c)**3) + 8*d/(3*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**2) - 2/((a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1507():
    f = 1/((a + b*x)**(sympy.S(7)/2)*(c + d*x)**(sympy.S(3)/2))
    F = -32*d**3*sqrt(a + b*x)/(5*sqrt(c + d*x)*(-a*d + b*c)**4) - 16*d**2/(5*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**3) + 4*d/(5*(a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-a*d + b*c)**2) - 2/((a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)*(-5*a*d + 5*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1508():
    f = 1/((a + b*x)**(sympy.S(9)/2)*(c + d*x)**(sympy.S(3)/2))
    F = 256*d**4*sqrt(a + b*x)/(35*sqrt(c + d*x)*(-a*d + b*c)**5) + 128*d**3/(35*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**4) - 32*d**2/(35*(a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-a*d + b*c)**3) + 16*d/(35*(a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)*(-a*d + b*c)**2) - 2/((a + b*x)**(sympy.S(7)/2)*sqrt(c + d*x)*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1509():
    f = 1/((a + b*x)**(sympy.S(11)/2)*(c + d*x)**(sympy.S(3)/2))
    F = -512*d**5*sqrt(a + b*x)/(63*sqrt(c + d*x)*(-a*d + b*c)**6) - 256*d**4/(63*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**5) + 64*d**3/(63*(a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-a*d + b*c)**4) - 32*d**2/(63*(a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)*(-a*d + b*c)**3) + 20*d/(63*(a + b*x)**(sympy.S(7)/2)*sqrt(c + d*x)*(-a*d + b*c)**2) - 2/((a + b*x)**(sympy.S(9)/2)*sqrt(c + d*x)*(-9*a*d + 9*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1510():
    f = (a + b*x)**(sympy.S(9)/2)/(c + d*x)**(sympy.S(5)/2)
    F = -105*b**(sympy.S(3)/2)*(-a*d + b*c)**3*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(8*d**(sympy.S(11)/2)) + 7*b**2*(a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)/d**3 - 35*b**2*(a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(-a*d + b*c)/(4*d**4) + 105*b**2*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)**2/(8*d**5) - 6*b*(a + b*x)**(sympy.S(7)/2)/(d**2*sqrt(c + d*x)) - 2*(a + b*x)**(sympy.S(9)/2)/(3*d*(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1511():
    f = (a + b*x)**(sympy.S(7)/2)/(c + d*x)**(sympy.S(5)/2)
    F = 35*b**(sympy.S(3)/2)*(-a*d + b*c)**2*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(4*d**(sympy.S(9)/2)) + 35*b**2*(a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)/(6*d**3) - 35*b**2*sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)/(4*d**4) - 14*b*(a + b*x)**(sympy.S(5)/2)/(3*d**2*sqrt(c + d*x)) - 2*(a + b*x)**(sympy.S(7)/2)/(3*d*(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1512():
    f = (a + b*x)**(sympy.S(5)/2)/(c + d*x)**(sympy.S(5)/2)
    F = -5*b**(sympy.S(3)/2)*(-a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/d**(sympy.S(7)/2) + 5*b**2*sqrt(a + b*x)*sqrt(c + d*x)/d**3 - 10*b*(a + b*x)**(sympy.S(3)/2)/(3*d**2*sqrt(c + d*x)) - 2*(a + b*x)**(sympy.S(5)/2)/(3*d*(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1513():
    f = (a + b*x)**(sympy.S(3)/2)/(c + d*x)**(sympy.S(5)/2)
    F = 2*b**(sympy.S(3)/2)*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/d**(sympy.S(5)/2) - 2*b*sqrt(a + b*x)/(d**2*sqrt(c + d*x)) - 2*(a + b*x)**(sympy.S(3)/2)/(3*d*(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1514():
    f = sqrt(a + b*x)/(c + d*x)**(sympy.S(5)/2)
    F = 2*(a + b*x)**(sympy.S(3)/2)/((c + d*x)**(sympy.S(3)/2)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1515():
    f = 1/(sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/2))
    F = 4*b*sqrt(a + b*x)/(3*sqrt(c + d*x)*(-a*d + b*c)**2) + 2*sqrt(a + b*x)/((c + d*x)**(sympy.S(3)/2)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1516():
    f = 1/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(5)/2))
    F = -16*b*d*sqrt(a + b*x)/(3*sqrt(c + d*x)*(-a*d + b*c)**3) - 8*d*sqrt(a + b*x)/(3*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**2) - 2/(sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1517():
    f = 1/((a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(5)/2))
    F = 32*b*d**2*sqrt(a + b*x)/(3*sqrt(c + d*x)*(-a*d + b*c)**4) + 16*d**2*sqrt(a + b*x)/(3*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**3) + 4*d/(sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**2) - 2/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1518():
    f = 1/((a + b*x)**(sympy.S(7)/2)*(c + d*x)**(sympy.S(5)/2))
    F = -256*b*d**3*sqrt(a + b*x)/(15*sqrt(c + d*x)*(-a*d + b*c)**5) - 128*d**3*sqrt(a + b*x)/(15*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**4) - 32*d**2/(5*sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**3) + 16*d/(15*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**2) - 2/((a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(3)/2)*(-5*a*d + 5*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1519():
    f = 1/((a + b*x)**(sympy.S(9)/2)*(c + d*x)**(sympy.S(5)/2))
    F = 512*b*d**4*sqrt(a + b*x)/(21*sqrt(c + d*x)*(-a*d + b*c)**6) + 256*d**4*sqrt(a + b*x)/(21*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**5) + 64*d**3/(7*sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**4) - 32*d**2/(21*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**3) + 4*d/(7*(a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**2) - 2/((a + b*x)**(sympy.S(7)/2)*(c + d*x)**(sympy.S(3)/2)*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1520():
    f = 1/(sqrt(a + b*x)*sqrt(a + b*x + 4))
    F = 2*asinh(sqrt(a + b*x)/2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1521():
    f = 1/(sqrt(b*x + 2)*sqrt(b*x + 6))
    F = 2*asinh(sqrt(b*x + 2)/2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1522():
    f = 1/(sqrt(b*x + 1)*sqrt(b*x + 5))
    F = 2*asinh(sqrt(b*x + 1)/2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1523():
    f = 1/(sqrt(b*x)*sqrt(b*x + 4))
    F = 2*asinh(sqrt(b*x)/2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1524():
    f = 1/(sqrt(b*x - 1)*sqrt(b*x + 3))
    F = 2*asinh(sqrt(b*x - 1)/2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1525():
    f = 1/(sqrt(b*x - 2)*sqrt(b*x + 2))
    F = acosh(b*x/2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1526():
    f = 1/(sqrt(b*x - 3)*sqrt(b*x + 1))
    F = 2*asinh(sqrt(b*x - 3)/2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1527():
    f = 1/(sqrt(b*x + 2)*sqrt(b*x + 3))
    F = 2*asinh(sqrt(b*x + 2))/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1528():
    f = 1/(b*x + 2)
    F = log(b*x + 2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1529():
    f = 1/(sqrt(b*x + 1)*sqrt(b*x + 2))
    F = 2*asinh(sqrt(b*x + 1))/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1530():
    f = 1/(sqrt(b*x)*sqrt(b*x + 2))
    F = 2*asinh(sqrt(2)*sqrt(b*x)/2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1531():
    f = 1/(sqrt(b*x - 1)*sqrt(b*x + 2))
    F = 2*asinh(sqrt(3)*sqrt(b*x - 1)/3)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1532():
    f = 1/(sqrt(b*x - 2)*sqrt(b*x + 2))
    F = acosh(b*x/2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1533():
    f = 1/(sqrt(b*x - 3)*sqrt(b*x + 2))
    F = 2*asinh(sqrt(5)*sqrt(b*x - 3)/5)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1534():
    f = 1/(sqrt(-b*x + 3)*sqrt(b*x + 2))
    F = asin(2*b*x/5 + sympy.S(-1)/5)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1535():
    f = 1/(sqrt(-b*x + 2)*sqrt(b*x + 2))
    F = asin(b*x/2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1536():
    f = 1/(sqrt(-b*x + 1)*sqrt(b*x + 2))
    F = asin(2*b*x/3 + sympy.S(1)/3)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1537():
    f = 1/(sqrt(-b*x)*sqrt(b*x + 2))
    F = asin(b*x + 1)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1538():
    f = 1/(sqrt(-b*x - 1)*sqrt(b*x + 2))
    F = asin(2*b*x + 3)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1539():
    f = 1/(sqrt(-b*x - 2)*sqrt(b*x + 2))
    F = sqrt(b*x + 2)*log(b*x + 2)/(b*sqrt(-b*x - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1540():
    f = 1/(sqrt(-b*x - 3)*sqrt(b*x + 2))
    F = -2*atan(sqrt(-b*x - 3)/sqrt(b*x + 2))/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1541():
    f = 1/(sqrt(-b*x + 2)*sqrt(-b*x + 3))
    F = -2*asinh(sqrt(-b*x + 2))/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1542():
    f = 1/(-b*x + 2)
    F = -log(-b*x + 2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1543():
    f = 1/(sqrt(-b*x + 1)*sqrt(-b*x + 2))
    F = -2*asinh(sqrt(-b*x + 1))/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1544():
    f = 1/(sqrt(-b*x)*sqrt(-b*x + 2))
    F = -2*asinh(sqrt(2)*sqrt(-b*x)/2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1545():
    f = 1/(sqrt(-b*x - 1)*sqrt(-b*x + 2))
    F = -2*asinh(sqrt(3)*sqrt(-b*x - 1)/3)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1546():
    f = 1/(sqrt(-b*x - 2)*sqrt(-b*x + 2))
    F = -acosh(-b*x/2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1547():
    f = 1/(sqrt(-b*x - 3)*sqrt(-b*x + 2))
    F = -2*asinh(sqrt(5)*sqrt(-b*x - 3)/5)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1548():
    f = 1/(sqrt(b*x - 4)*sqrt(b*x + 4))
    F = acosh(b*x/4)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1549():
    f = 1/(sqrt(c + d*x)*sqrt(b*x + (b*c - b)/d))
    F = 2*asinh(sqrt(d)*sqrt(b*x - b*(1 - c)/d)/sqrt(b))/(sqrt(b)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1550():
    f = 1/(sqrt(x)*sqrt(2*x - 3))
    F = sqrt(2)*asinh(sqrt(3)*sqrt(2*x - 3)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1551():
    f = 1/(sqrt(2*x - 3)*sqrt(3*x + 2))
    F = sqrt(6)*asinh(sqrt(39)*sqrt(2*x - 3)/13)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1552():
    f = 1/(sqrt(c - d*x)*sqrt(b*x + (-b*c + b)/d))
    F = 2*asin(sqrt(d)*sqrt(b*x + b*(1 - c)/d)/sqrt(b))/(sqrt(b)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1553():
    f = 1/(sqrt(x)*sqrt(4 - x))
    F = asin(x/2 - 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1554():
    f = 1/(sqrt(x)*sqrt(3 - 2*x))
    F = sqrt(2)*asin(sqrt(6)*sqrt(x)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1555():
    f = 1/(sqrt(3 - 2*x)*sqrt(5*x + 3))
    F = sqrt(10)*asin(sqrt(42)*sqrt(5*x + 3)/21)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1556():
    f = 1/(sqrt(a - b*x)*sqrt(c + d*x))
    F = -2*atan(sqrt(d)*sqrt(a - b*x)/(sqrt(b)*sqrt(c + d*x)))/(sqrt(b)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1557():
    f = (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/3)
    F = 6*(a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(1)/3)/(17*b) + (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/3)*(-12*a*d + 12*b*c)/(187*b*d) - 108*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**2/(935*b*d**2) - 108*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(-a*d + b*c)**3*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(935*b**(sympy.S(4)/3)*d**3*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1558():
    f = sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/3)
    F = 6*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/3)/(11*b) + sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/3)*(-12*a*d + 12*b*c)/(55*b*d) + 12*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(-a*d + b*c)**2*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(55*b**(sympy.S(4)/3)*d**2*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1559():
    f = (c + d*x)**(sympy.S(1)/3)/sqrt(a + b*x)
    F = 6*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/3)/(5*b) - 4*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(-a*d + b*c)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(5*b**(sympy.S(4)/3)*d*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1560():
    f = (c + d*x)**(sympy.S(1)/3)/(a + b*x)**(sympy.S(3)/2)
    F = -2*(c + d*x)**(sympy.S(1)/3)/(b*sqrt(a + b*x)) - 4*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(3*b**(sympy.S(4)/3)*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1561():
    f = (c + d*x)**(sympy.S(1)/3)/(a + b*x)**(sympy.S(5)/2)
    F = -4*d*(c + d*x)**(sympy.S(1)/3)/(9*b*sqrt(a + b*x)*(-a*d + b*c)) - 2*(c + d*x)**(sympy.S(1)/3)/(3*b*(a + b*x)**(sympy.S(3)/2)) + 4*3**(sympy.S(3)/4)*d*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(27*b**(sympy.S(4)/3)*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1562():
    f = (c + d*x)**(sympy.S(1)/3)/(a + b*x)**(sympy.S(7)/2)
    F = 28*d**2*(c + d*x)**(sympy.S(1)/3)/(135*b*sqrt(a + b*x)*(-a*d + b*c)**2) - 4*d*(c + d*x)**(sympy.S(1)/3)/(45*b*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)) - 2*(c + d*x)**(sympy.S(1)/3)/(5*b*(a + b*x)**(sympy.S(5)/2)) - 28*3**(sympy.S(3)/4)*d**2*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(405*b**(sympy.S(4)/3)*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1563():
    f = (a + b*x)**(sympy.S(3)/2)/(c + d*x)**(sympy.S(1)/3)
    F = 6*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(2)/3)/(13*d) - sqrt(a + b*x)*(c + d*x)**(sympy.S(2)/3)*(-54*a*d + 54*b*c)/(91*d**2) - 162*sqrt(a + b*x)*(-a*d + b*c)**2/(91*b**(sympy.S(2)/3)*d**2*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))) + 81*3**(sympy.S(1)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(-a*d + b*c)**(sympy.S(7)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(91*b**(sympy.S(2)/3)*d**3*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) - 54*sqrt(2)*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*(-a*d + b*c)**(sympy.S(7)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(91*b**(sympy.S(2)/3)*d**3*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1564():
    f = sqrt(a + b*x)/(c + d*x)**(sympy.S(1)/3)
    F = 6*sqrt(a + b*x)*(c + d*x)**(sympy.S(2)/3)/(7*d) + sqrt(a + b*x)*(-18*a*d + 18*b*c)/(7*b**(sympy.S(2)/3)*d*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))) - 9*3**(sympy.S(1)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(-a*d + b*c)**(sympy.S(4)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(7*b**(sympy.S(2)/3)*d**2*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) + 6*sqrt(2)*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*(-a*d + b*c)**(sympy.S(4)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(7*b**(sympy.S(2)/3)*d**2*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1565():
    f = 1/(sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/3))
    F = -6*sqrt(a + b*x)/(b**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))) + 3*3**(sympy.S(1)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(b**(sympy.S(2)/3)*d*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) - 2*sqrt(2)*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(b**(sympy.S(2)/3)*d*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1566():
    f = 1/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/3))
    F = -2*(c + d*x)**(sympy.S(2)/3)/(sqrt(a + b*x)*(-a*d + b*c)) - 2*d*sqrt(a + b*x)/(b**(sympy.S(2)/3)*(-a*d + b*c)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))) + 3**(sympy.S(1)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(2)/3)) - 2*sqrt(2)*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1567():
    f = 1/((a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(1)/3))
    F = 10*d*(c + d*x)**(sympy.S(2)/3)/(9*sqrt(a + b*x)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(2)/3)/((a + b*x)**(sympy.S(3)/2)*(-3*a*d + 3*b*c)) + 10*d**2*sqrt(a + b*x)/(9*b**(sympy.S(2)/3)*(-a*d + b*c)**2*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))) - 5*3**(sympy.S(1)/4)*d*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(9*b**(sympy.S(2)/3)*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/3)) + 10*sqrt(2)*3**(sympy.S(3)/4)*d*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(27*b**(sympy.S(2)/3)*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1568():
    f = (a + b*x)**(sympy.S(3)/2)/(c + d*x)**(sympy.S(2)/3)
    F = 6*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/3)/(11*d) - sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/3)*(-54*a*d + 54*b*c)/(55*d**2) - 54*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(-a*d + b*c)**2*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(55*b**(sympy.S(1)/3)*d**3*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1569():
    f = sqrt(a + b*x)/(c + d*x)**(sympy.S(2)/3)
    F = 6*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/3)/(5*d) + 6*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(-a*d + b*c)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(5*b**(sympy.S(1)/3)*d**2*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1570():
    f = 1/(sqrt(a + b*x)*(c + d*x)**(sympy.S(2)/3))
    F = -2*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(b**(sympy.S(1)/3)*d*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1571():
    f = 1/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(2)/3))
    F = -2*(c + d*x)**(sympy.S(1)/3)/(sqrt(a + b*x)*(-a*d + b*c)) + 2*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(3*b**(sympy.S(1)/3)*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1572():
    f = 1/((a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(2)/3))
    F = 14*d*(c + d*x)**(sympy.S(1)/3)/(9*sqrt(a + b*x)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(1)/3)/((a + b*x)**(sympy.S(3)/2)*(-3*a*d + 3*b*c)) - 14*3**(sympy.S(3)/4)*d*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(asin((-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(27*b**(sympy.S(1)/3)*sqrt(-(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1573():
    f = (a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(1)/3)
    F = (a + b*x)**(sympy.S(5)/3)*(c + d*x)**(sympy.S(1)/3)/(2*b) + (a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)/(6*b*d) + (-a*d + b*c)**2*log(-1 + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)))/(6*b**(sympy.S(4)/3)*d**(sympy.S(5)/3)) + (-a*d + b*c)**2*log(c + d*x)/(18*b**(sympy.S(4)/3)*d**(sympy.S(5)/3)) + sqrt(3)*(-a*d + b*c)**2*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(3*b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)))/(9*b**(sympy.S(4)/3)*d**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1574():
    f = (c + d*x)**(sympy.S(1)/3)/(a + b*x)**(sympy.S(1)/3)
    F = (a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(1)/3)/b - (-a*d + b*c)*log(-1 + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)))/(2*b**(sympy.S(4)/3)*d**(sympy.S(2)/3)) - (-a*d + b*c)*log(c + d*x)/(6*b**(sympy.S(4)/3)*d**(sympy.S(2)/3)) - sqrt(3)*(-a*d + b*c)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(3*b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)))/(3*b**(sympy.S(4)/3)*d**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1575():
    f = (c + d*x)**(sympy.S(1)/3)/(a + b*x)**(sympy.S(4)/3)
    F = -3*(c + d*x)**(sympy.S(1)/3)/(b*(a + b*x)**(sympy.S(1)/3)) - 3*d**(sympy.S(1)/3)*log(-1 + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)))/(2*b**(sympy.S(4)/3)) - d**(sympy.S(1)/3)*log(c + d*x)/(2*b**(sympy.S(4)/3)) - sqrt(3)*d**(sympy.S(1)/3)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(3*b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)))/b**(sympy.S(4)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1576():
    f = (c + d*x)**(sympy.S(1)/3)/(a + b*x)**(sympy.S(7)/3)
    F = -3*(c + d*x)**(sympy.S(4)/3)/((a + b*x)**(sympy.S(4)/3)*(-4*a*d + 4*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1577():
    f = (c + d*x)**(sympy.S(1)/3)/(a + b*x)**(sympy.S(10)/3)
    F = 9*d*(c + d*x)**(sympy.S(4)/3)/(28*(a + b*x)**(sympy.S(4)/3)*(-a*d + b*c)**2) - 3*(c + d*x)**(sympy.S(4)/3)/((a + b*x)**(sympy.S(7)/3)*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1578():
    f = (c + d*x)**(sympy.S(1)/3)/(a + b*x)**(sympy.S(13)/3)
    F = -27*d**2*(c + d*x)**(sympy.S(4)/3)/(140*(a + b*x)**(sympy.S(4)/3)*(-a*d + b*c)**3) + 9*d*(c + d*x)**(sympy.S(4)/3)/(35*(a + b*x)**(sympy.S(7)/3)*(-a*d + b*c)**2) - 3*(c + d*x)**(sympy.S(4)/3)/((a + b*x)**(sympy.S(10)/3)*(-10*a*d + 10*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1579():
    f = (c + d*x)**(sympy.S(1)/3)/(a + b*x)**(sympy.S(16)/3)
    F = 243*d**3*(c + d*x)**(sympy.S(4)/3)/(1820*(a + b*x)**(sympy.S(4)/3)*(-a*d + b*c)**4) - 81*d**2*(c + d*x)**(sympy.S(4)/3)/(455*(a + b*x)**(sympy.S(7)/3)*(-a*d + b*c)**3) + 27*d*(c + d*x)**(sympy.S(4)/3)/(130*(a + b*x)**(sympy.S(10)/3)*(-a*d + b*c)**2) - 3*(c + d*x)**(sympy.S(4)/3)/((a + b*x)**(sympy.S(13)/3)*(-13*a*d + 13*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1580():
    f = (a + b*x)**(sympy.S(4)/3)*(c + d*x)**(sympy.S(1)/3)
    F = 3*(a + b*x)**(sympy.S(7)/3)*(c + d*x)**(sympy.S(1)/3)/(8*b) + (a + b*x)**(sympy.S(4)/3)*(c + d*x)**(sympy.S(1)/3)*(-3*a*d + 3*b*c)/(40*b*d) - 3*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**2/(20*b*d**2) + 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*((a + b*x)*(c + d*x))**(sympy.S(2)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(sqrt(3) + 2)*(-a*d + b*c)**3*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(20*b**(sympy.S(4)/3)*d**(sympy.S(7)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1581():
    f = (a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)
    F = 3*(a + b*x)**(sympy.S(4)/3)*(c + d*x)**(sympy.S(1)/3)/(5*b) + (a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-3*a*d + 3*b*c)/(10*b*d) - 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*((a + b*x)*(c + d*x))**(sympy.S(2)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(sqrt(3) + 2)*(-a*d + b*c)**2*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(10*b**(sympy.S(4)/3)*d**(sympy.S(4)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1582():
    f = (c + d*x)**(sympy.S(1)/3)/(a + b*x)**(sympy.S(2)/3)
    F = 3*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(2*b) + 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*((a + b*x)*(c + d*x))**(sympy.S(2)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(sqrt(3) + 2)*(-a*d + b*c)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(2*b**(sympy.S(4)/3)*d**(sympy.S(1)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1583():
    f = (c + d*x)**(sympy.S(1)/3)/(a + b*x)**(sympy.S(5)/3)
    F = -3*(c + d*x)**(sympy.S(1)/3)/(2*b*(a + b*x)**(sympy.S(2)/3)) + 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(sqrt(3) + 2)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(2*b**(sympy.S(4)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1584():
    f = (c + d*x)**(sympy.S(1)/3)/(a + b*x)**(sympy.S(8)/3)
    F = -3*d*(c + d*x)**(sympy.S(1)/3)/(10*b*(a + b*x)**(sympy.S(2)/3)*(-a*d + b*c)) - 3*(c + d*x)**(sympy.S(1)/3)/(5*b*(a + b*x)**(sympy.S(5)/3)) - 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*d**(sympy.S(5)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(sqrt(3) + 2)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(10*b**(sympy.S(4)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)*(-a*d + b*c)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1585():
    f = (a + b*x)**(sympy.S(4)/3)/(c + d*x)**(sympy.S(1)/3)
    F = (a + b*x)**(sympy.S(4)/3)*(c + d*x)**(sympy.S(2)/3)/(2*d) + (a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(2)/3)*(2*a*d - 2*b*c)/(3*d**2) - (-a*d + b*c)**2*log(a + b*x)/(9*b**(sympy.S(2)/3)*d**(sympy.S(7)/3)) - (-a*d + b*c)**2*log(b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)) - 1)/(3*b**(sympy.S(2)/3)*d**(sympy.S(7)/3)) - 2*sqrt(3)*(-a*d + b*c)**2*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(3*d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)) + sqrt(3)/3)/(9*b**(sympy.S(2)/3)*d**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1586():
    f = (a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3)
    F = (a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(2)/3)/d + (-a*d + b*c)*log(a + b*x)/(6*b**(sympy.S(2)/3)*d**(sympy.S(4)/3)) + (-a*d + b*c)*log(b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)) - 1)/(2*b**(sympy.S(2)/3)*d**(sympy.S(4)/3)) + sqrt(3)*(-a*d + b*c)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(3*d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)) + sqrt(3)/3)/(3*b**(sympy.S(2)/3)*d**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1587():
    f = 1/((a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(1)/3))
    F = -log(a + b*x)/(2*b**(sympy.S(2)/3)*d**(sympy.S(1)/3)) - 3*log(b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)) - 1)/(2*b**(sympy.S(2)/3)*d**(sympy.S(1)/3)) - sqrt(3)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(3*d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)) + sqrt(3)/3)/(b**(sympy.S(2)/3)*d**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1588():
    f = 1/((a + b*x)**(sympy.S(5)/3)*(c + d*x)**(sympy.S(1)/3))
    F = -3*(c + d*x)**(sympy.S(2)/3)/((a + b*x)**(sympy.S(2)/3)*(-2*a*d + 2*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1589():
    f = 1/((a + b*x)**(sympy.S(8)/3)*(c + d*x)**(sympy.S(1)/3))
    F = 9*d*(c + d*x)**(sympy.S(2)/3)/(10*(a + b*x)**(sympy.S(2)/3)*(-a*d + b*c)**2) - 3*(c + d*x)**(sympy.S(2)/3)/((a + b*x)**(sympy.S(5)/3)*(-5*a*d + 5*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1590():
    f = 1/((a + b*x)**(sympy.S(11)/3)*(c + d*x)**(sympy.S(1)/3))
    F = -27*d**2*(c + d*x)**(sympy.S(2)/3)/(40*(a + b*x)**(sympy.S(2)/3)*(-a*d + b*c)**3) + 9*d*(c + d*x)**(sympy.S(2)/3)/(20*(a + b*x)**(sympy.S(5)/3)*(-a*d + b*c)**2) - 3*(c + d*x)**(sympy.S(2)/3)/((a + b*x)**(sympy.S(8)/3)*(-8*a*d + 8*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1591():
    f = 1/((a + b*x)**(sympy.S(14)/3)*(c + d*x)**(sympy.S(1)/3))
    F = 243*d**3*(c + d*x)**(sympy.S(2)/3)/(440*(a + b*x)**(sympy.S(2)/3)*(-a*d + b*c)**4) - 81*d**2*(c + d*x)**(sympy.S(2)/3)/(220*(a + b*x)**(sympy.S(5)/3)*(-a*d + b*c)**3) + 27*d*(c + d*x)**(sympy.S(2)/3)/(88*(a + b*x)**(sympy.S(8)/3)*(-a*d + b*c)**2) - 3*(c + d*x)**(sympy.S(2)/3)/((a + b*x)**(sympy.S(11)/3)*(-11*a*d + 11*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1592():
    f = (a + b*x)**(sympy.S(8)/3)/(c + d*x)**(sympy.S(1)/3)
    F = 3*(a + b*x)**(sympy.S(8)/3)*(c + d*x)**(sympy.S(2)/3)/(10*d) - (a + b*x)**(sympy.S(5)/3)*(c + d*x)**(sympy.S(2)/3)*(-12*a*d + 12*b*c)/(35*d**2) + 3*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)*(-a*d + b*c)**2/(7*d**3) + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(2 - sqrt(3))*(-a*d + b*c)**(sympy.S(11)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(14*b**(sympy.S(2)/3)*d**(sympy.S(11)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 2*2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(-a*d + b*c)**(sympy.S(11)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(7*b**(sympy.S(2)/3)*d**(sympy.S(11)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 3*2**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**3*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(7*b**(sympy.S(2)/3)*d**(sympy.S(11)/3)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))*(a*d + b*c + 2*b*d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1593():
    f = (a + b*x)**(sympy.S(5)/3)/(c + d*x)**(sympy.S(1)/3)
    F = 3*(a + b*x)**(sympy.S(5)/3)*(c + d*x)**(sympy.S(2)/3)/(7*d) - (a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)*(-15*a*d + 15*b*c)/(28*d**2) - 15*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(2 - sqrt(3))*(-a*d + b*c)**(sympy.S(8)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(56*b**(sympy.S(2)/3)*d**(sympy.S(8)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 5*2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(-a*d + b*c)**(sympy.S(8)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(14*b**(sympy.S(2)/3)*d**(sympy.S(8)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 15*2**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**2*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(28*b**(sympy.S(2)/3)*d**(sympy.S(8)/3)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))*(a*d + b*c + 2*b*d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1594():
    f = (a + b*x)**(sympy.S(2)/3)/(c + d*x)**(sympy.S(1)/3)
    F = 3*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)/(4*d) + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(2 - sqrt(3))*(-a*d + b*c)**(sympy.S(5)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(8*b**(sympy.S(2)/3)*d**(sympy.S(5)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(-a*d + b*c)**(sympy.S(5)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(2*b**(sympy.S(2)/3)*d**(sympy.S(5)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 2**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-3*a*d + 3*b*c)*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(4*b**(sympy.S(2)/3)*d**(sympy.S(5)/3)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))*(a*d + b*c + 2*b*d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1595():
    f = 1/((a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3))
    F = -3*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(2 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(4*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 3*2**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(2*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))*(a*d + b*c + 2*b*d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1596():
    f = 1/((a + b*x)**(sympy.S(4)/3)*(c + d*x)**(sympy.S(1)/3))
    F = -3*(c + d*x)**(sympy.S(2)/3)/((a + b*x)**(sympy.S(1)/3)*(-a*d + b*c)) - 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(2 - sqrt(3))*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(4*b**(sympy.S(2)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 3*2**(sympy.S(2)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(2*b**(sympy.S(2)/3)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))*(a*d + b*c + 2*b*d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1597():
    f = 1/((a + b*x)**(sympy.S(7)/3)*(c + d*x)**(sympy.S(1)/3))
    F = 3*d*(c + d*x)**(sympy.S(2)/3)/(2*(a + b*x)**(sympy.S(1)/3)*(-a*d + b*c)**2) - 3*(c + d*x)**(sympy.S(2)/3)/((a + b*x)**(sympy.S(4)/3)*(-4*a*d + 4*b*c)) + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*d**(sympy.S(4)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(2 - sqrt(3))*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(8*b**(sympy.S(2)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(4)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*d**(sympy.S(4)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(2*b**(sympy.S(2)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(4)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 3*2**(sympy.S(2)/3)*d**(sympy.S(4)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(4*b**(sympy.S(2)/3)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**2*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))*(a*d + b*c + 2*b*d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1598():
    f = 1/((a + b*x)**(sympy.S(10)/3)*(c + d*x)**(sympy.S(1)/3))
    F = -15*d**2*(c + d*x)**(sympy.S(2)/3)/(14*(a + b*x)**(sympy.S(1)/3)*(-a*d + b*c)**3) + 15*d*(c + d*x)**(sympy.S(2)/3)/(28*(a + b*x)**(sympy.S(4)/3)*(-a*d + b*c)**2) - 3*(c + d*x)**(sympy.S(2)/3)/((a + b*x)**(sympy.S(7)/3)*(-7*a*d + 7*b*c)) - 15*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*d**(sympy.S(7)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(2 - sqrt(3))*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(56*b**(sympy.S(2)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(7)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 5*2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*d**(sympy.S(7)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(14*b**(sympy.S(2)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(7)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 15*2**(sympy.S(2)/3)*d**(sympy.S(7)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(28*b**(sympy.S(2)/3)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**3*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))*(a*d + b*c + 2*b*d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1599():
    f = (a + b*x)**(sympy.S(5)/3)/(c + d*x)**(sympy.S(2)/3)
    F = (a + b*x)**(sympy.S(5)/3)*(c + d*x)**(sympy.S(1)/3)/(2*d) + (a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(1)/3)*(5*a*d - 5*b*c)/(6*d**2) - 5*(-a*d + b*c)**2*log(-1 + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)))/(6*b**(sympy.S(1)/3)*d**(sympy.S(8)/3)) - 5*(-a*d + b*c)**2*log(c + d*x)/(18*b**(sympy.S(1)/3)*d**(sympy.S(8)/3)) - 5*sqrt(3)*(-a*d + b*c)**2*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(3*b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)))/(9*b**(sympy.S(1)/3)*d**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1600():
    f = (a + b*x)**(sympy.S(2)/3)/(c + d*x)**(sympy.S(2)/3)
    F = (a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(1)/3)/d + sqrt(3)*(-2*a*d + 2*b*c)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(3*b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)))/(3*b**(sympy.S(1)/3)*d**(sympy.S(5)/3)) + (-a*d + b*c)*log(-1 + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)))/(b**(sympy.S(1)/3)*d**(sympy.S(5)/3)) + (-a*d + b*c)*log(c + d*x)/(3*b**(sympy.S(1)/3)*d**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1601():
    f = 1/((a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(2)/3))
    F = -3*log(-1 + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)))/(2*b**(sympy.S(1)/3)*d**(sympy.S(2)/3)) - log(c + d*x)/(2*b**(sympy.S(1)/3)*d**(sympy.S(2)/3)) - sqrt(3)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(3*b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)))/(b**(sympy.S(1)/3)*d**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1602():
    f = 1/((a + b*x)**(sympy.S(4)/3)*(c + d*x)**(sympy.S(2)/3))
    F = -3*(c + d*x)**(sympy.S(1)/3)/((a + b*x)**(sympy.S(1)/3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1603():
    f = 1/((a + b*x)**(sympy.S(7)/3)*(c + d*x)**(sympy.S(2)/3))
    F = 9*d*(c + d*x)**(sympy.S(1)/3)/(4*(a + b*x)**(sympy.S(1)/3)*(-a*d + b*c)**2) - 3*(c + d*x)**(sympy.S(1)/3)/((a + b*x)**(sympy.S(4)/3)*(-4*a*d + 4*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1604():
    f = 1/((a + b*x)**(sympy.S(10)/3)*(c + d*x)**(sympy.S(2)/3))
    F = -27*d**2*(c + d*x)**(sympy.S(1)/3)/(14*(a + b*x)**(sympy.S(1)/3)*(-a*d + b*c)**3) + 9*d*(c + d*x)**(sympy.S(1)/3)/(14*(a + b*x)**(sympy.S(4)/3)*(-a*d + b*c)**2) - 3*(c + d*x)**(sympy.S(1)/3)/((a + b*x)**(sympy.S(7)/3)*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1605():
    f = 1/((a + b*x)**(sympy.S(13)/3)*(c + d*x)**(sympy.S(2)/3))
    F = 243*d**3*(c + d*x)**(sympy.S(1)/3)/(140*(a + b*x)**(sympy.S(1)/3)*(-a*d + b*c)**4) - 81*d**2*(c + d*x)**(sympy.S(1)/3)/(140*(a + b*x)**(sympy.S(4)/3)*(-a*d + b*c)**3) + 27*d*(c + d*x)**(sympy.S(1)/3)/(70*(a + b*x)**(sympy.S(7)/3)*(-a*d + b*c)**2) - 3*(c + d*x)**(sympy.S(1)/3)/((a + b*x)**(sympy.S(10)/3)*(-10*a*d + 10*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1606():
    f = (a + b*x)**(sympy.S(7)/3)/(c + d*x)**(sympy.S(2)/3)
    F = 3*(a + b*x)**(sympy.S(7)/3)*(c + d*x)**(sympy.S(1)/3)/(8*d) - (a + b*x)**(sympy.S(4)/3)*(c + d*x)**(sympy.S(1)/3)*(-21*a*d + 21*b*c)/(40*d**2) + 21*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**2/(20*d**3) - 7*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*((a + b*x)*(c + d*x))**(sympy.S(2)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(sqrt(3) + 2)*(-a*d + b*c)**3*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(20*b**(sympy.S(1)/3)*d**(sympy.S(10)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1607():
    f = (a + b*x)**(sympy.S(4)/3)/(c + d*x)**(sympy.S(2)/3)
    F = 3*(a + b*x)**(sympy.S(4)/3)*(c + d*x)**(sympy.S(1)/3)/(5*d) - (a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-6*a*d + 6*b*c)/(5*d**2) + 2*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*((a + b*x)*(c + d*x))**(sympy.S(2)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(sqrt(3) + 2)*(-a*d + b*c)**2*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(5*b**(sympy.S(1)/3)*d**(sympy.S(7)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1608():
    f = (a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(2)/3)
    F = 3*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(2*d) - 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*((a + b*x)*(c + d*x))**(sympy.S(2)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(sqrt(3) + 2)*(-a*d + b*c)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(2*b**(sympy.S(1)/3)*d**(sympy.S(4)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1609():
    f = 1/((a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3))
    F = 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*((a + b*x)*(c + d*x))**(sympy.S(2)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(sqrt(3) + 2)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1610():
    f = 1/((a + b*x)**(sympy.S(5)/3)*(c + d*x)**(sympy.S(2)/3))
    F = -3*(c + d*x)**(sympy.S(1)/3)/((a + b*x)**(sympy.S(2)/3)*(-2*a*d + 2*b*c)) - 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(sqrt(3) + 2)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(2*b**(sympy.S(1)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)*(-a*d + b*c)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1611():
    f = 1/((a + b*x)**(sympy.S(8)/3)*(c + d*x)**(sympy.S(2)/3))
    F = 6*d*(c + d*x)**(sympy.S(1)/3)/(5*(a + b*x)**(sympy.S(2)/3)*(-a*d + b*c)**2) - 3*(c + d*x)**(sympy.S(1)/3)/((a + b*x)**(sympy.S(5)/3)*(-5*a*d + 5*b*c)) + 2*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*d**(sympy.S(5)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(sqrt(3) + 2)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(5*b**(sympy.S(1)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)*(-a*d + b*c)**2*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1612():
    f = 1/((a + b*x)**(sympy.S(11)/3)*(c + d*x)**(sympy.S(2)/3))
    F = -21*d**2*(c + d*x)**(sympy.S(1)/3)/(20*(a + b*x)**(sympy.S(2)/3)*(-a*d + b*c)**3) + 21*d*(c + d*x)**(sympy.S(1)/3)/(40*(a + b*x)**(sympy.S(5)/3)*(-a*d + b*c)**2) - 3*(c + d*x)**(sympy.S(1)/3)/((a + b*x)**(sympy.S(8)/3)*(-8*a*d + 8*b*c)) - 7*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*d**(sympy.S(8)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(sqrt(3) + 2)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(20*b**(sympy.S(1)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)*(-a*d + b*c)**3*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1613():
    f = (a + b*x)**(sympy.S(7)/3)/(c + d*x)**(sympy.S(4)/3)
    F = -7*b**(sympy.S(1)/3)*(-a*d + b*c)**2*log(a + b*x)/(9*d**(sympy.S(10)/3)) - 7*b**(sympy.S(1)/3)*(-a*d + b*c)**2*log(b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)) - 1)/(3*d**(sympy.S(10)/3)) - 14*sqrt(3)*b**(sympy.S(1)/3)*(-a*d + b*c)**2*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(3*d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)) + sqrt(3)/3)/(9*d**(sympy.S(10)/3)) + 7*b*(a + b*x)**(sympy.S(4)/3)*(c + d*x)**(sympy.S(2)/3)/(2*d**2) - 14*b*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(2)/3)*(-a*d + b*c)/(3*d**3) - 3*(a + b*x)**(sympy.S(7)/3)/(d*(c + d*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1614():
    f = (a + b*x)**(sympy.S(4)/3)/(c + d*x)**(sympy.S(4)/3)
    F = 2*b**(sympy.S(1)/3)*(-a*d + b*c)*log(a + b*x)/(3*d**(sympy.S(7)/3)) + 2*b**(sympy.S(1)/3)*(-a*d + b*c)*log(b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)) - 1)/d**(sympy.S(7)/3) + 4*sqrt(3)*b**(sympy.S(1)/3)*(-a*d + b*c)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(3*d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)) + sqrt(3)/3)/(3*d**(sympy.S(7)/3)) + 4*b*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(2)/3)/d**2 - 3*(a + b*x)**(sympy.S(4)/3)/(d*(c + d*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1615():
    f = (a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(4)/3)
    F = -b**(sympy.S(1)/3)*log(a + b*x)/(2*d**(sympy.S(4)/3)) - 3*b**(sympy.S(1)/3)*log(b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)) - 1)/(2*d**(sympy.S(4)/3)) - sqrt(3)*b**(sympy.S(1)/3)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)/(3*d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)) + sqrt(3)/3)/d**(sympy.S(4)/3) - 3*(a + b*x)**(sympy.S(1)/3)/(d*(c + d*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1616():
    f = 1/((a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(4)/3))
    F = 3*(a + b*x)**(sympy.S(1)/3)/((c + d*x)**(sympy.S(1)/3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1617():
    f = 1/((a + b*x)**(sympy.S(5)/3)*(c + d*x)**(sympy.S(4)/3))
    F = -9*d*(a + b*x)**(sympy.S(1)/3)/(2*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**2) - 3/((a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(1)/3)*(-2*a*d + 2*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1618():
    f = 1/((a + b*x)**(sympy.S(8)/3)*(c + d*x)**(sympy.S(4)/3))
    F = 27*d**2*(a + b*x)**(sympy.S(1)/3)/(5*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**3) + 9*d/(5*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**2) - 3/((a + b*x)**(sympy.S(5)/3)*(c + d*x)**(sympy.S(1)/3)*(-5*a*d + 5*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1619():
    f = 1/((a + b*x)**(sympy.S(11)/3)*(c + d*x)**(sympy.S(4)/3))
    F = -243*d**3*(a + b*x)**(sympy.S(1)/3)/(40*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**4) - 81*d**2/(40*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**3) + 27*d/(40*(a + b*x)**(sympy.S(5)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**2) - 3/((a + b*x)**(sympy.S(8)/3)*(c + d*x)**(sympy.S(1)/3)*(-8*a*d + 8*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1620():
    f = (a + b*x)**(sympy.S(8)/3)/(c + d*x)**(sympy.S(4)/3)
    F = -15*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(2 - sqrt(3))*(-a*d + b*c)**(sympy.S(8)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(7*d**(sympy.S(11)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 20*2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(-a*d + b*c)**(sympy.S(8)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(7*d**(sympy.S(11)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 30*2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**2*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(7*d**(sympy.S(11)/3)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))*(a*d + b*c + 2*b*d*x)) + 24*b*(a + b*x)**(sympy.S(5)/3)*(c + d*x)**(sympy.S(2)/3)/(7*d**2) - 30*b*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)*(-a*d + b*c)/(7*d**3) - 3*(a + b*x)**(sympy.S(8)/3)/(d*(c + d*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1621():
    f = (a + b*x)**(sympy.S(5)/3)/(c + d*x)**(sympy.S(4)/3)
    F = 15*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(2 - sqrt(3))*(-a*d + b*c)**(sympy.S(5)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(8*d**(sympy.S(8)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 5*2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(-a*d + b*c)**(sympy.S(5)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(2*d**(sympy.S(8)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 15*2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(4*d**(sympy.S(8)/3)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))*(a*d + b*c + 2*b*d*x)) + 15*b*(a + b*x)**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3)/(4*d**2) - 3*(a + b*x)**(sympy.S(5)/3)/(d*(c + d*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1622():
    f = (a + b*x)**(sympy.S(2)/3)/(c + d*x)**(sympy.S(4)/3)
    F = -3*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(2 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(2*d**(sympy.S(5)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 2*2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(d**(sympy.S(5)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 3*2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(d**(sympy.S(5)/3)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))*(a*d + b*c + 2*b*d*x)) - 3*(a + b*x)**(sympy.S(2)/3)/(d*(c + d*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1623():
    f = 1/((a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(4)/3))
    F = 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(2 - sqrt(3))*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(4*d**(sympy.S(2)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(d**(sympy.S(2)/3)*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 3*2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(2*d**(sympy.S(2)/3)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))*(a*d + b*c + 2*b*d*x)) + 3*(a + b*x)**(sympy.S(2)/3)/((c + d*x)**(sympy.S(1)/3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1624():
    f = 1/((a + b*x)**(sympy.S(4)/3)*(c + d*x)**(sympy.S(4)/3))
    F = -3*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(2 - sqrt(3))*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(2*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(4)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 2*2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(4)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 3*2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/((a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**2*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))*(a*d + b*c + 2*b*d*x)) - 6*d*(a + b*x)**(sympy.S(2)/3)/((c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**2) - 3/((a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1625():
    f = 1/((a + b*x)**(sympy.S(7)/3)*(c + d*x)**(sympy.S(4)/3))
    F = 15*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*d**(sympy.S(4)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*sqrt(2 - sqrt(3))*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(8*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(7)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 5*2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*d**(sympy.S(4)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((2*2**(sympy.S(1)/3)*b**(sympy.S(2)/3)*d**(sympy.S(2)/3)*((a + b*x)*(c + d*x))**(sympy.S(2)/3) - 2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(2)/3) + (-a*d + b*c)**(sympy.S(4)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(asin((2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 - sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))), -7 - 4*sqrt(3))/(2*sqrt((-a*d + b*c)**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))**2)*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(7)/3)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 15*2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(4)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3)*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(4*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**3*(2**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d**(sympy.S(1)/3)*((a + b*x)*(c + d*x))**(sympy.S(1)/3) + (1 + sqrt(3))*(-a*d + b*c)**(sympy.S(2)/3))*(a*d + b*c + 2*b*d*x)) + 15*d**2*(a + b*x)**(sympy.S(2)/3)/(2*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**3) + 15*d/(4*(a + b*x)**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**2) - 3/((a + b*x)**(sympy.S(4)/3)*(c + d*x)**(sympy.S(1)/3)*(-4*a*d + 4*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1626():
    f = (x - 1)**(sympy.S(1)/3)/(x + 1)**(sympy.S(1)/3)
    F = (x - 1)**(sympy.S(1)/3)*(x + 1)**(sympy.S(2)/3) + log(-1 + (x + 1)**(sympy.S(1)/3)/(x - 1)**(sympy.S(1)/3)) + log(x - 1)/3 + 2*sqrt(3)*atan(sqrt(3)/3 + 2*sqrt(3)*(x + 1)**(sympy.S(1)/3)/(3*(x - 1)**(sympy.S(1)/3)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1627():
    f = (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/4)
    F = 4*(a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(1)/4)/(11*b) + (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/4)*(-4*a*d + 4*b*c)/(77*b*d) - 8*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2/(77*b*d**2) + 16*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(13)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(77*b**(sympy.S(5)/4)*d**3*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1628():
    f = sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/4)
    F = 4*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/4)/(7*b) + sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/4)*(-4*a*d + 4*b*c)/(21*b*d) - 8*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(9)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(21*b**(sympy.S(5)/4)*d**2*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1629():
    f = (c + d*x)**(sympy.S(1)/4)/sqrt(a + b*x)
    F = 4*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/4)/(3*b) + 4*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(5)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(3*b**(sympy.S(5)/4)*d*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1630():
    f = (c + d*x)**(sympy.S(1)/4)/(a + b*x)**(sympy.S(3)/2)
    F = -2*(c + d*x)**(sympy.S(1)/4)/(b*sqrt(a + b*x)) + 2*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(1)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(b**(sympy.S(5)/4)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1631():
    f = (c + d*x)**(sympy.S(1)/4)/(a + b*x)**(sympy.S(5)/2)
    F = -d*(c + d*x)**(sympy.S(1)/4)/(3*b*sqrt(a + b*x)*(-a*d + b*c)) - 2*(c + d*x)**(sympy.S(1)/4)/(3*b*(a + b*x)**(sympy.S(3)/2)) - d*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(3*b**(sympy.S(5)/4)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1632():
    f = (c + d*x)**(sympy.S(1)/4)/(a + b*x)**(sympy.S(7)/2)
    F = d**2*(c + d*x)**(sympy.S(1)/4)/(6*b*sqrt(a + b*x)*(-a*d + b*c)**2) - d*(c + d*x)**(sympy.S(1)/4)/(15*b*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)) - 2*(c + d*x)**(sympy.S(1)/4)/(5*b*(a + b*x)**(sympy.S(5)/2)) + d**2*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(6*b**(sympy.S(5)/4)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1633():
    f = (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/4)
    F = 4*(a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(3)/4)/(13*b) + (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/4)*(-4*a*d + 4*b*c)/(39*b*d) - 8*sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/4)*(-a*d + b*c)**2/(65*b*d**2) + 16*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(15)/4)*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(65*b**(sympy.S(7)/4)*d**3*sqrt(a + b*x)) - 16*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(15)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(65*b**(sympy.S(7)/4)*d**3*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1634():
    f = sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/4)
    F = 4*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/4)/(9*b) + sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/4)*(-4*a*d + 4*b*c)/(15*b*d) - 8*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(11)/4)*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(15*b**(sympy.S(7)/4)*d**2*sqrt(a + b*x)) + 8*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(11)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(15*b**(sympy.S(7)/4)*d**2*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1635():
    f = (c + d*x)**(sympy.S(3)/4)/sqrt(a + b*x)
    F = 4*sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/4)/(5*b) + 12*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(7)/4)*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*b**(sympy.S(7)/4)*d*sqrt(a + b*x)) - 12*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(7)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*b**(sympy.S(7)/4)*d*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1636():
    f = (c + d*x)**(sympy.S(3)/4)/(a + b*x)**(sympy.S(3)/2)
    F = -2*(c + d*x)**(sympy.S(3)/4)/(b*sqrt(a + b*x)) + 6*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(3)/4)*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(b**(sympy.S(7)/4)*sqrt(a + b*x)) - 6*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(3)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(b**(sympy.S(7)/4)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1637():
    f = (c + d*x)**(sympy.S(3)/4)/(a + b*x)**(sympy.S(5)/2)
    F = -d*(c + d*x)**(sympy.S(3)/4)/(b*sqrt(a + b*x)*(-a*d + b*c)) - 2*(c + d*x)**(sympy.S(3)/4)/(3*b*(a + b*x)**(sympy.S(3)/2)) + d*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(b**(sympy.S(7)/4)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(1)/4)) - d*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(b**(sympy.S(7)/4)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1638():
    f = (c + d*x)**(sympy.S(3)/4)/(a + b*x)**(sympy.S(7)/2)
    F = 3*d**2*(c + d*x)**(sympy.S(3)/4)/(10*b*sqrt(a + b*x)*(-a*d + b*c)**2) - d*(c + d*x)**(sympy.S(3)/4)/(5*b*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)) - 2*(c + d*x)**(sympy.S(3)/4)/(5*b*(a + b*x)**(sympy.S(5)/2)) - 3*d**2*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(10*b**(sympy.S(7)/4)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/4)) + 3*d**2*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(10*b**(sympy.S(7)/4)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1639():
    f = (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(5)/4)
    F = 4*(a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(5)/4)/(15*b) + (a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(1)/4)*(-4*a*d + 4*b*c)/(33*b**2) + 4*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2/(231*b**2*d) - 8*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**3/(231*b**2*d**2) + 16*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(17)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(231*b**(sympy.S(9)/4)*d**3*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1640():
    f = sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/4)
    F = 4*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(5)/4)/(11*b) + (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/4)*(-20*a*d + 20*b*c)/(77*b**2) + 20*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2/(231*b**2*d) - 40*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(13)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(231*b**(sympy.S(9)/4)*d**2*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1641():
    f = (c + d*x)**(sympy.S(5)/4)/sqrt(a + b*x)
    F = 4*sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/4)/(7*b) + sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/4)*(-20*a*d + 20*b*c)/(21*b**2) + 20*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(9)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(21*b**(sympy.S(9)/4)*d*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1642():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(3)/2)
    F = -2*(c + d*x)**(sympy.S(5)/4)/(b*sqrt(a + b*x)) + 10*d*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/4)/(3*b**2) + 10*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(5)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(3*b**(sympy.S(9)/4)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1643():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(5)/2)
    F = -2*(c + d*x)**(sympy.S(5)/4)/(3*b*(a + b*x)**(sympy.S(3)/2)) - 5*d*(c + d*x)**(sympy.S(1)/4)/(3*b**2*sqrt(a + b*x)) + 5*d*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(1)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(3*b**(sympy.S(9)/4)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1644():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(7)/2)
    F = -2*(c + d*x)**(sympy.S(5)/4)/(5*b*(a + b*x)**(sympy.S(5)/2)) - d**2*(c + d*x)**(sympy.S(1)/4)/(6*b**2*sqrt(a + b*x)*(-a*d + b*c)) - d*(c + d*x)**(sympy.S(1)/4)/(3*b**2*(a + b*x)**(sympy.S(3)/2)) - d**2*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(6*b**(sympy.S(9)/4)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1645():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(9)/2)
    F = -2*(c + d*x)**(sympy.S(5)/4)/(7*b*(a + b*x)**(sympy.S(7)/2)) + 5*d**3*(c + d*x)**(sympy.S(1)/4)/(84*b**2*sqrt(a + b*x)*(-a*d + b*c)**2) - d**2*(c + d*x)**(sympy.S(1)/4)/(42*b**2*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)) - d*(c + d*x)**(sympy.S(1)/4)/(7*b**2*(a + b*x)**(sympy.S(5)/2)) + 5*d**3*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(84*b**(sympy.S(9)/4)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1646():
    f = (a + b*x)**(sympy.S(5)/2)/(c + d*x)**(sympy.S(1)/4)
    F = 4*(a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(3)/4)/(13*d) - (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/4)*(-40*a*d + 40*b*c)/(117*d**2) + 16*sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/4)*(-a*d + b*c)**2/(39*d**3) - 32*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(15)/4)*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(39*b**(sympy.S(3)/4)*d**4*sqrt(a + b*x)) + 32*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(15)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(39*b**(sympy.S(3)/4)*d**4*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1647():
    f = (a + b*x)**(sympy.S(3)/2)/(c + d*x)**(sympy.S(1)/4)
    F = 4*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/4)/(9*d) + sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/4)*(8*a*d - 8*b*c)/(15*d**2) + 16*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(11)/4)*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(15*b**(sympy.S(3)/4)*d**3*sqrt(a + b*x)) - 16*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(11)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(15*b**(sympy.S(3)/4)*d**3*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1648():
    f = sqrt(a + b*x)/(c + d*x)**(sympy.S(1)/4)
    F = 4*sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/4)/(5*d) - 8*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(7)/4)*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*b**(sympy.S(3)/4)*d**2*sqrt(a + b*x)) + 8*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(7)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*b**(sympy.S(3)/4)*d**2*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1649():
    f = 1/(sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/4))
    F = 4*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(3)/4)*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(b**(sympy.S(3)/4)*d*sqrt(a + b*x)) - 4*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(3)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(b**(sympy.S(3)/4)*d*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1650():
    f = 1/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/4))
    F = -2*(c + d*x)**(sympy.S(3)/4)/(sqrt(a + b*x)*(-a*d + b*c)) + 2*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(b**(sympy.S(3)/4)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(1)/4)) - 2*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(b**(sympy.S(3)/4)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1651():
    f = 1/((a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(1)/4))
    F = d*(c + d*x)**(sympy.S(3)/4)/(sqrt(a + b*x)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(3)/4)/((a + b*x)**(sympy.S(3)/2)*(-3*a*d + 3*b*c)) - d*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(b**(sympy.S(3)/4)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/4)) + d*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(b**(sympy.S(3)/4)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1652():
    f = (a + b*x)**(sympy.S(3)/2)/(c + d*x)**(sympy.S(3)/4)
    F = 4*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/4)/(7*d) + sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/4)*(8*a*d - 8*b*c)/(7*d**2) + 16*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(9)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(7*b**(sympy.S(1)/4)*d**3*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1653():
    f = sqrt(a + b*x)/(c + d*x)**(sympy.S(3)/4)
    F = 4*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/4)/(3*d) - 8*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(5)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(3*b**(sympy.S(1)/4)*d**2*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1654():
    f = 1/(sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/4))
    F = 4*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(1)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(b**(sympy.S(1)/4)*d*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1655():
    f = 1/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/4))
    F = -2*(c + d*x)**(sympy.S(1)/4)/(sqrt(a + b*x)*(-a*d + b*c)) - 2*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(b**(sympy.S(1)/4)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1656():
    f = 1/((a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(3)/4))
    F = 5*d*(c + d*x)**(sympy.S(1)/4)/(3*sqrt(a + b*x)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(1)/4)/((a + b*x)**(sympy.S(3)/2)*(-3*a*d + 3*b*c)) + 5*d*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(3*b**(sympy.S(1)/4)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1657():
    f = (a + b*x)**(sympy.S(5)/2)/(c + d*x)**(sympy.S(5)/4)
    F = 32*b**(sympy.S(1)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(11)/4)*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(3*d**4*sqrt(a + b*x)) - 32*b**(sympy.S(1)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(11)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(3*d**4*sqrt(a + b*x)) + 40*b*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/4)/(9*d**2) - 16*b*sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/4)*(-a*d + b*c)/(3*d**3) - 4*(a + b*x)**(sympy.S(5)/2)/(d*(c + d*x)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1658():
    f = (a + b*x)**(sympy.S(3)/2)/(c + d*x)**(sympy.S(5)/4)
    F = -48*b**(sympy.S(1)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(7)/4)*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*d**3*sqrt(a + b*x)) + 48*b**(sympy.S(1)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(7)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*d**3*sqrt(a + b*x)) + 24*b*sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/4)/(5*d**2) - 4*(a + b*x)**(sympy.S(3)/2)/(d*(c + d*x)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1659():
    f = sqrt(a + b*x)/(c + d*x)**(sympy.S(5)/4)
    F = 8*b**(sympy.S(1)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(3)/4)*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(d**2*sqrt(a + b*x)) - 8*b**(sympy.S(1)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(3)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(d**2*sqrt(a + b*x)) - 4*sqrt(a + b*x)/(d*(c + d*x)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1660():
    f = 1/(sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/4))
    F = -4*b**(sympy.S(1)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(d*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(1)/4)) + 4*b**(sympy.S(1)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(d*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(1)/4)) + 4*sqrt(a + b*x)/((c + d*x)**(sympy.S(1)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1661():
    f = 1/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(5)/4))
    F = 6*b**(sympy.S(1)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/4)) - 6*b**(sympy.S(1)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/4)) - 6*d*sqrt(a + b*x)/((c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2) - 2/(sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1662():
    f = 1/((a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(5)/4))
    F = -7*b**(sympy.S(1)/4)*d*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(9)/4)) + 7*b**(sympy.S(1)/4)*d*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(9)/4)) + 7*d**2*sqrt(a + b*x)/((c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**3) + 7*d/(3*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2) - 2/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/4)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1663():
    f = (a + b*x)**(sympy.S(7)/2)/(c + d*x)**(sympy.S(7)/4)
    F = -320*b**(sympy.S(3)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(13)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(33*d**5*sqrt(a + b*x)) + 56*b*(a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(1)/4)/(33*d**2) - 80*b*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)/(33*d**3) + 160*b*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2/(33*d**4) - 4*(a + b*x)**(sympy.S(7)/2)/(3*d*(c + d*x)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1664():
    f = (a + b*x)**(sympy.S(3)/2)/(c + d*x)**(sympy.S(7)/4)
    F = -16*b**(sympy.S(3)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(5)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(3*d**3*sqrt(a + b*x)) + 8*b*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/4)/(3*d**2) - 4*(a + b*x)**(sympy.S(3)/2)/(3*d*(c + d*x)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1665():
    f = sqrt(a + b*x)/(c + d*x)**(sympy.S(7)/4)
    F = 8*b**(sympy.S(3)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(1)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(3*d**2*sqrt(a + b*x)) - 4*sqrt(a + b*x)/(3*d*(c + d*x)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1666():
    f = 1/(sqrt(a + b*x)*(c + d*x)**(sympy.S(7)/4))
    F = 4*b**(sympy.S(3)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(3*d*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(3)/4)) + 4*sqrt(a + b*x)/((c + d*x)**(sympy.S(3)/4)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1667():
    f = 1/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(7)/4))
    F = -10*b**(sympy.S(3)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(3*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(7)/4)) - 10*d*sqrt(a + b*x)/(3*(c + d*x)**(sympy.S(3)/4)*(-a*d + b*c)**2) - 2/(sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1668():
    f = 1/((a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(7)/4))
    F = 5*b**(sympy.S(3)/4)*d*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(11)/4)) + 5*d**2*sqrt(a + b*x)/((c + d*x)**(sympy.S(3)/4)*(-a*d + b*c)**3) + 3*d/(sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/4)*(-a*d + b*c)**2) - 2/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/4)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1669():
    f = (a + b*x)**(sympy.S(7)/2)/(c + d*x)**(sympy.S(9)/4)
    F = 448*b**(sympy.S(5)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(11)/4)*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(15*d**5*sqrt(a + b*x)) - 448*b**(sympy.S(5)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(11)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(15*d**5*sqrt(a + b*x)) + 112*b**2*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/4)/(9*d**3) - 224*b**2*sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/4)*(-a*d + b*c)/(15*d**4) - 56*b*(a + b*x)**(sympy.S(5)/2)/(5*d**2*(c + d*x)**(sympy.S(1)/4)) - 4*(a + b*x)**(sympy.S(7)/2)/(5*d*(c + d*x)**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1670():
    f = (a + b*x)**(sympy.S(5)/2)/(c + d*x)**(sympy.S(9)/4)
    F = -96*b**(sympy.S(5)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(7)/4)*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*d**4*sqrt(a + b*x)) + 96*b**(sympy.S(5)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(7)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*d**4*sqrt(a + b*x)) + 48*b**2*sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/4)/(5*d**3) - 8*b*(a + b*x)**(sympy.S(3)/2)/(d**2*(c + d*x)**(sympy.S(1)/4)) - 4*(a + b*x)**(sympy.S(5)/2)/(5*d*(c + d*x)**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1671():
    f = (a + b*x)**(sympy.S(3)/2)/(c + d*x)**(sympy.S(9)/4)
    F = 48*b**(sympy.S(5)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(3)/4)*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*d**3*sqrt(a + b*x)) - 48*b**(sympy.S(5)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*(-a*d + b*c)**(sympy.S(3)/4)*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*d**3*sqrt(a + b*x)) - 24*b*sqrt(a + b*x)/(5*d**2*(c + d*x)**(sympy.S(1)/4)) - 4*(a + b*x)**(sympy.S(3)/2)/(5*d*(c + d*x)**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1672():
    f = sqrt(a + b*x)/(c + d*x)**(sympy.S(9)/4)
    F = -8*b**(sympy.S(5)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*d**2*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(1)/4)) + 8*b**(sympy.S(5)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*d**2*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(1)/4)) + 8*b*sqrt(a + b*x)/(5*d*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)) - 4*sqrt(a + b*x)/(5*d*(c + d*x)**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1673():
    f = 1/(sqrt(a + b*x)*(c + d*x)**(sympy.S(9)/4))
    F = -12*b**(sympy.S(5)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*d*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/4)) + 12*b**(sympy.S(5)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*d*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/4)) + 12*b*sqrt(a + b*x)/(5*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2) + 4*sqrt(a + b*x)/((c + d*x)**(sympy.S(5)/4)*(-5*a*d + 5*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1674():
    f = 1/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(9)/4))
    F = 42*b**(sympy.S(5)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(9)/4)) - 42*b**(sympy.S(5)/4)*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(9)/4)) - 42*b*d*sqrt(a + b*x)/(5*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**3) - 14*d*sqrt(a + b*x)/(5*(c + d*x)**(sympy.S(5)/4)*(-a*d + b*c)**2) - 2/(sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1675():
    f = 1/((a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(9)/4))
    F = -77*b**(sympy.S(5)/4)*d*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_e(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(13)/4)) + 77*b**(sympy.S(5)/4)*d*sqrt(-d*(a + b*x)/(-a*d + b*c))*elliptic_f(asin(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(-a*d + b*c)**(sympy.S(1)/4)), -1)/(5*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(13)/4)) + 77*b*d**2*sqrt(a + b*x)/(5*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**4) + 77*d**2*sqrt(a + b*x)/(15*(c + d*x)**(sympy.S(5)/4)*(-a*d + b*c)**3) + 11*d/(3*sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/4)*(-a*d + b*c)**2) - 2/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(5)/4)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1676():
    f = (a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(5)/4)
    F = (a + b*x)**(sympy.S(7)/4)*(c + d*x)**(sympy.S(5)/4)/(3*b) + (a + b*x)**(sympy.S(7)/4)*(c + d*x)**(sympy.S(1)/4)*(-5*a*d + 5*b*c)/(24*b**2) + 5*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2/(96*b**2*d) + 5*(-a*d + b*c)**3*atan(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(64*b**(sympy.S(9)/4)*d**(sympy.S(7)/4)) - 5*(-a*d + b*c)**3*atanh(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(64*b**(sympy.S(9)/4)*d**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1677():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(1)/4)
    F = (a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(5)/4)/(2*b) + (a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(1)/4)*(-5*a*d + 5*b*c)/(8*b**2) - 5*(-a*d + b*c)**2*atan(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(16*b**(sympy.S(9)/4)*d**(sympy.S(3)/4)) + 5*(-a*d + b*c)**2*atanh(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(16*b**(sympy.S(9)/4)*d**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1678():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(5)/4)
    F = -4*(c + d*x)**(sympy.S(5)/4)/(b*(a + b*x)**(sympy.S(1)/4)) + 5*d*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(1)/4)/b**2 - 5*d**(sympy.S(1)/4)*(-a*d + b*c)*atan(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(2*b**(sympy.S(9)/4)) + 5*d**(sympy.S(1)/4)*(-a*d + b*c)*atanh(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(2*b**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1679():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(9)/4)
    F = -4*(c + d*x)**(sympy.S(5)/4)/(5*b*(a + b*x)**(sympy.S(5)/4)) - 4*d*(c + d*x)**(sympy.S(1)/4)/(b**2*(a + b*x)**(sympy.S(1)/4)) - 2*d**(sympy.S(5)/4)*atan(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/b**(sympy.S(9)/4) + 2*d**(sympy.S(5)/4)*atanh(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/b**(sympy.S(9)/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1680():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(13)/4)
    F = -4*(c + d*x)**(sympy.S(9)/4)/((a + b*x)**(sympy.S(9)/4)*(-9*a*d + 9*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1681():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(17)/4)
    F = 16*d*(c + d*x)**(sympy.S(9)/4)/(117*(a + b*x)**(sympy.S(9)/4)*(-a*d + b*c)**2) - 4*(c + d*x)**(sympy.S(9)/4)/((a + b*x)**(sympy.S(13)/4)*(-13*a*d + 13*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1682():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(21)/4)
    F = -128*d**2*(c + d*x)**(sympy.S(9)/4)/(1989*(a + b*x)**(sympy.S(9)/4)*(-a*d + b*c)**3) + 32*d*(c + d*x)**(sympy.S(9)/4)/(221*(a + b*x)**(sympy.S(13)/4)*(-a*d + b*c)**2) - 4*(c + d*x)**(sympy.S(9)/4)/((a + b*x)**(sympy.S(17)/4)*(-17*a*d + 17*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1683():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(25)/4)
    F = 512*d**3*(c + d*x)**(sympy.S(9)/4)/(13923*(a + b*x)**(sympy.S(9)/4)*(-a*d + b*c)**4) - 128*d**2*(c + d*x)**(sympy.S(9)/4)/(1547*(a + b*x)**(sympy.S(13)/4)*(-a*d + b*c)**3) + 16*d*(c + d*x)**(sympy.S(9)/4)/(119*(a + b*x)**(sympy.S(17)/4)*(-a*d + b*c)**2) - 4*(c + d*x)**(sympy.S(9)/4)/((a + b*x)**(sympy.S(21)/4)*(-21*a*d + 21*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1684():
    f = (a + b*x)**(sympy.S(5)/4)*(c + d*x)**(sympy.S(5)/4)
    F = 2*(a + b*x)**(sympy.S(9)/4)*(c + d*x)**(sympy.S(5)/4)/(7*b) + (a + b*x)**(sympy.S(9)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)/(7*b**2) + (a + b*x)**(sympy.S(5)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2/(42*b**2*d) - 5*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**3/(84*b**2*d**2) + 5*sqrt(2)*((a + b*x)*(c + d*x))**(sympy.S(3)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(9)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(336*b**(sympy.S(9)/4)*d**(sympy.S(9)/4)*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1685():
    f = (a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(5)/4)
    F = 2*(a + b*x)**(sympy.S(5)/4)*(c + d*x)**(sympy.S(5)/4)/(5*b) + (a + b*x)**(sympy.S(5)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)/(3*b**2) + (a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2/(6*b**2*d) - sqrt(2)*((a + b*x)*(c + d*x))**(sympy.S(3)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(7)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(24*b**(sympy.S(9)/4)*d**(sympy.S(5)/4)*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1686():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(3)/4)
    F = 2*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(5)/4)/(3*b) + (a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(-5*a*d + 5*b*c)/(3*b**2) + 5*sqrt(2)*((a + b*x)*(c + d*x))**(sympy.S(3)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(5)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(12*b**(sympy.S(9)/4)*d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1687():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(7)/4)
    F = -4*(c + d*x)**(sympy.S(5)/4)/(3*b*(a + b*x)**(sympy.S(3)/4)) + 10*d*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/(3*b**2) + 5*sqrt(2)*d**(sympy.S(3)/4)*((a + b*x)*(c + d*x))**(sympy.S(3)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(3)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(6*b**(sympy.S(9)/4)*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1688():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(11)/4)
    F = -4*(c + d*x)**(sympy.S(5)/4)/(7*b*(a + b*x)**(sympy.S(7)/4)) - 20*d*(c + d*x)**(sympy.S(1)/4)/(21*b**2*(a + b*x)**(sympy.S(3)/4)) + 5*sqrt(2)*d**(sympy.S(7)/4)*((a + b*x)*(c + d*x))**(sympy.S(3)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*sqrt(-a*d + b*c)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(21*b**(sympy.S(9)/4)*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1689():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(15)/4)
    F = -4*(c + d*x)**(sympy.S(5)/4)/(11*b*(a + b*x)**(sympy.S(11)/4)) - 20*d**2*(c + d*x)**(sympy.S(1)/4)/(231*b**2*(a + b*x)**(sympy.S(3)/4)*(-a*d + b*c)) - 20*d*(c + d*x)**(sympy.S(1)/4)/(77*b**2*(a + b*x)**(sympy.S(7)/4)) - 10*sqrt(2)*d**(sympy.S(11)/4)*((a + b*x)*(c + d*x))**(sympy.S(3)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(231*b**(sympy.S(9)/4)*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)*sqrt(-a*d + b*c)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1690():
    f = (c + d*x)**(sympy.S(5)/4)/(a + b*x)**(sympy.S(19)/4)
    F = -4*(c + d*x)**(sympy.S(5)/4)/(15*b*(a + b*x)**(sympy.S(15)/4)) + 8*d**3*(c + d*x)**(sympy.S(1)/4)/(231*b**2*(a + b*x)**(sympy.S(3)/4)*(-a*d + b*c)**2) - 4*d**2*(c + d*x)**(sympy.S(1)/4)/(231*b**2*(a + b*x)**(sympy.S(7)/4)*(-a*d + b*c)) - 4*d*(c + d*x)**(sympy.S(1)/4)/(33*b**2*(a + b*x)**(sympy.S(11)/4)) + 4*sqrt(2)*d**(sympy.S(15)/4)*((a + b*x)*(c + d*x))**(sympy.S(3)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(231*b**(sympy.S(9)/4)*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)*(-a*d + b*c)**(sympy.S(3)/2)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1691():
    f = (a + b*x)**(sympy.S(5)/4)/(c + d*x)**(sympy.S(1)/4)
    F = (a + b*x)**(sympy.S(5)/4)*(c + d*x)**(sympy.S(3)/4)/(2*d) + (a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(3)/4)*(5*a*d - 5*b*c)/(8*d**2) + 5*(-a*d + b*c)**2*atan(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(16*b**(sympy.S(3)/4)*d**(sympy.S(9)/4)) + 5*(-a*d + b*c)**2*atanh(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(16*b**(sympy.S(3)/4)*d**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1692():
    f = (a + b*x)**(sympy.S(1)/4)/(c + d*x)**(sympy.S(1)/4)
    F = (a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(3)/4)/d - (-a*d + b*c)*atan(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(2*b**(sympy.S(3)/4)*d**(sympy.S(5)/4)) - (-a*d + b*c)*atanh(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(2*b**(sympy.S(3)/4)*d**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1693():
    f = 1/((a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(1)/4))
    F = 2*atan(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(b**(sympy.S(3)/4)*d**(sympy.S(1)/4)) + 2*atanh(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(b**(sympy.S(3)/4)*d**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1694():
    f = 1/((a + b*x)**(sympy.S(7)/4)*(c + d*x)**(sympy.S(1)/4))
    F = -4*(c + d*x)**(sympy.S(3)/4)/((a + b*x)**(sympy.S(3)/4)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1695():
    f = 1/((a + b*x)**(sympy.S(11)/4)*(c + d*x)**(sympy.S(1)/4))
    F = 16*d*(c + d*x)**(sympy.S(3)/4)/(21*(a + b*x)**(sympy.S(3)/4)*(-a*d + b*c)**2) - 4*(c + d*x)**(sympy.S(3)/4)/((a + b*x)**(sympy.S(7)/4)*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1696():
    f = 1/((a + b*x)**(sympy.S(15)/4)*(c + d*x)**(sympy.S(1)/4))
    F = -128*d**2*(c + d*x)**(sympy.S(3)/4)/(231*(a + b*x)**(sympy.S(3)/4)*(-a*d + b*c)**3) + 32*d*(c + d*x)**(sympy.S(3)/4)/(77*(a + b*x)**(sympy.S(7)/4)*(-a*d + b*c)**2) - 4*(c + d*x)**(sympy.S(3)/4)/((a + b*x)**(sympy.S(11)/4)*(-11*a*d + 11*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1697():
    f = 1/((a + b*x)**(sympy.S(19)/4)*(c + d*x)**(sympy.S(1)/4))
    F = 512*d**3*(c + d*x)**(sympy.S(3)/4)/(1155*(a + b*x)**(sympy.S(3)/4)*(-a*d + b*c)**4) - 128*d**2*(c + d*x)**(sympy.S(3)/4)/(385*(a + b*x)**(sympy.S(7)/4)*(-a*d + b*c)**3) + 16*d*(c + d*x)**(sympy.S(3)/4)/(55*(a + b*x)**(sympy.S(11)/4)*(-a*d + b*c)**2) - 4*(c + d*x)**(sympy.S(3)/4)/((a + b*x)**(sympy.S(15)/4)*(-15*a*d + 15*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1698():
    f = (a + b*x)**(sympy.S(7)/4)/(c + d*x)**(sympy.S(1)/4)
    F = 2*(a + b*x)**(sympy.S(7)/4)*(c + d*x)**(sympy.S(3)/4)/(5*d) - (a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)*(-7*a*d + 7*b*c)/(15*d**2) + sqrt((a + b*x)*(c + d*x))*(-7*a*d + 7*b*c)*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(10*sqrt(b)*d**(sympy.S(5)/2)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*(a*d + b*c + 2*b*d*x)) - 7*sqrt(2)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(7)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(20*b**(sympy.S(3)/4)*d**(sympy.S(11)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 7*sqrt(2)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(7)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(40*b**(sympy.S(3)/4)*d**(sympy.S(11)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1699():
    f = (a + b*x)**(sympy.S(3)/4)/(c + d*x)**(sympy.S(1)/4)
    F = 2*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)/(3*d) - sqrt((a + b*x)*(c + d*x))*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(sqrt(b)*d**(sympy.S(3)/2)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*(a*d + b*c + 2*b*d*x)) + sqrt(2)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(5)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(2*b**(sympy.S(3)/4)*d**(sympy.S(7)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - sqrt(2)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(5)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(4*b**(sympy.S(3)/4)*d**(sympy.S(7)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1700():
    f = 1/((a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4))
    F = 2*sqrt((a + b*x)*(c + d*x))*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(sqrt(b)*sqrt(d)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*(a*d + b*c + 2*b*d*x)) - sqrt(2)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(3)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(b**(sympy.S(3)/4)*d**(sympy.S(3)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + sqrt(2)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(3)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(2*b**(sympy.S(3)/4)*d**(sympy.S(3)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1701():
    f = 1/((a + b*x)**(sympy.S(5)/4)*(c + d*x)**(sympy.S(1)/4))
    F = -4*(c + d*x)**(sympy.S(3)/4)/((a + b*x)**(sympy.S(1)/4)*(-a*d + b*c)) + 4*sqrt(d)*sqrt((a + b*x)*(c + d*x))*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(sqrt(b)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*(a*d + b*c + 2*b*d*x)) - 2*sqrt(2)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*sqrt(-a*d + b*c)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(b**(sympy.S(3)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + sqrt(2)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*sqrt(-a*d + b*c)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(b**(sympy.S(3)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1702():
    f = 1/((a + b*x)**(sympy.S(9)/4)*(c + d*x)**(sympy.S(1)/4))
    F = 8*d*(c + d*x)**(sympy.S(3)/4)/(5*(a + b*x)**(sympy.S(1)/4)*(-a*d + b*c)**2) - 4*(c + d*x)**(sympy.S(3)/4)/((a + b*x)**(sympy.S(5)/4)*(-5*a*d + 5*b*c)) - 8*d**(sympy.S(3)/2)*sqrt((a + b*x)*(c + d*x))*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(5*sqrt(b)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**3*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*(a*d + b*c + 2*b*d*x)) + 4*sqrt(2)*d**(sympy.S(5)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*sqrt(-a*d + b*c)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 2*sqrt(2)*d**(sympy.S(5)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*sqrt(-a*d + b*c)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1703():
    f = (a + b*x)**(sympy.S(7)/4)/(c + d*x)**(sympy.S(3)/4)
    F = (a + b*x)**(sympy.S(7)/4)*(c + d*x)**(sympy.S(1)/4)/(2*d) + (a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(1)/4)*(7*a*d - 7*b*c)/(8*d**2) - 21*(-a*d + b*c)**2*atan(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(16*b**(sympy.S(1)/4)*d**(sympy.S(11)/4)) + 21*(-a*d + b*c)**2*atanh(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(16*b**(sympy.S(1)/4)*d**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1704():
    f = (a + b*x)**(sympy.S(3)/4)/(c + d*x)**(sympy.S(3)/4)
    F = (a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(1)/4)/d + (-3*a*d + 3*b*c)*atan(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(2*b**(sympy.S(1)/4)*d**(sympy.S(7)/4)) - (-3*a*d + 3*b*c)*atanh(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(2*b**(sympy.S(1)/4)*d**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1705():
    f = 1/((a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(3)/4))
    F = -2*atan(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(b**(sympy.S(1)/4)*d**(sympy.S(3)/4)) + 2*atanh(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(b**(sympy.S(1)/4)*d**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1706():
    f = 1/((a + b*x)**(sympy.S(5)/4)*(c + d*x)**(sympy.S(3)/4))
    F = -4*(c + d*x)**(sympy.S(1)/4)/((a + b*x)**(sympy.S(1)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1707():
    f = 1/((a + b*x)**(sympy.S(9)/4)*(c + d*x)**(sympy.S(3)/4))
    F = 16*d*(c + d*x)**(sympy.S(1)/4)/(5*(a + b*x)**(sympy.S(1)/4)*(-a*d + b*c)**2) - 4*(c + d*x)**(sympy.S(1)/4)/((a + b*x)**(sympy.S(5)/4)*(-5*a*d + 5*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1708():
    f = 1/((a + b*x)**(sympy.S(13)/4)*(c + d*x)**(sympy.S(3)/4))
    F = -128*d**2*(c + d*x)**(sympy.S(1)/4)/(45*(a + b*x)**(sympy.S(1)/4)*(-a*d + b*c)**3) + 32*d*(c + d*x)**(sympy.S(1)/4)/(45*(a + b*x)**(sympy.S(5)/4)*(-a*d + b*c)**2) - 4*(c + d*x)**(sympy.S(1)/4)/((a + b*x)**(sympy.S(9)/4)*(-9*a*d + 9*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1709():
    f = 1/((a + b*x)**(sympy.S(17)/4)*(c + d*x)**(sympy.S(3)/4))
    F = 512*d**3*(c + d*x)**(sympy.S(1)/4)/(195*(a + b*x)**(sympy.S(1)/4)*(-a*d + b*c)**4) - 128*d**2*(c + d*x)**(sympy.S(1)/4)/(195*(a + b*x)**(sympy.S(5)/4)*(-a*d + b*c)**3) + 16*d*(c + d*x)**(sympy.S(1)/4)/(39*(a + b*x)**(sympy.S(9)/4)*(-a*d + b*c)**2) - 4*(c + d*x)**(sympy.S(1)/4)/((a + b*x)**(sympy.S(13)/4)*(-13*a*d + 13*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1710():
    f = (a + b*x)**(sympy.S(5)/4)/(c + d*x)**(sympy.S(3)/4)
    F = 2*(a + b*x)**(sympy.S(5)/4)*(c + d*x)**(sympy.S(1)/4)/(3*d) - (a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(-5*a*d + 5*b*c)/(3*d**2) + 5*sqrt(2)*((a + b*x)*(c + d*x))**(sympy.S(3)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(5)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(12*b**(sympy.S(1)/4)*d**(sympy.S(9)/4)*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1711():
    f = (a + b*x)**(sympy.S(1)/4)/(c + d*x)**(sympy.S(3)/4)
    F = 2*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)/d - sqrt(2)*((a + b*x)*(c + d*x))**(sympy.S(3)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(3)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(2*b**(sympy.S(1)/4)*d**(sympy.S(5)/4)*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1712():
    f = 1/((a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4))
    F = sqrt(2)*((a + b*x)*(c + d*x))**(sympy.S(3)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*sqrt(-a*d + b*c)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1713():
    f = 1/((a + b*x)**(sympy.S(7)/4)*(c + d*x)**(sympy.S(3)/4))
    F = -4*(c + d*x)**(sympy.S(1)/4)/((a + b*x)**(sympy.S(3)/4)*(-3*a*d + 3*b*c)) - 2*sqrt(2)*d**(sympy.S(3)/4)*((a + b*x)*(c + d*x))**(sympy.S(3)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(3*b**(sympy.S(1)/4)*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)*sqrt(-a*d + b*c)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1714():
    f = 1/((a + b*x)**(sympy.S(11)/4)*(c + d*x)**(sympy.S(3)/4))
    F = 8*d*(c + d*x)**(sympy.S(1)/4)/(7*(a + b*x)**(sympy.S(3)/4)*(-a*d + b*c)**2) - 4*(c + d*x)**(sympy.S(1)/4)/((a + b*x)**(sympy.S(7)/4)*(-7*a*d + 7*b*c)) + 4*sqrt(2)*d**(sympy.S(7)/4)*((a + b*x)*(c + d*x))**(sympy.S(3)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(7*b**(sympy.S(1)/4)*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)*(-a*d + b*c)**(sympy.S(3)/2)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1715():
    f = (a + b*x)**(sympy.S(5)/4)/(c + d*x)**(sympy.S(5)/4)
    F = -5*b**(sympy.S(1)/4)*(-a*d + b*c)*atan(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(2*d**(sympy.S(9)/4)) - 5*b**(sympy.S(1)/4)*(-a*d + b*c)*atanh(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/(2*d**(sympy.S(9)/4)) + 5*b*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(3)/4)/d**2 - 4*(a + b*x)**(sympy.S(5)/4)/(d*(c + d*x)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1716():
    f = (a + b*x)**(sympy.S(1)/4)/(c + d*x)**(sympy.S(5)/4)
    F = 2*b**(sympy.S(1)/4)*atan(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/d**(sympy.S(5)/4) + 2*b**(sympy.S(1)/4)*atanh(d**(sympy.S(1)/4)*(a + b*x)**(sympy.S(1)/4)/(b**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)))/d**(sympy.S(5)/4) - 4*(a + b*x)**(sympy.S(1)/4)/(d*(c + d*x)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1717():
    f = 1/((a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(5)/4))
    F = 4*(a + b*x)**(sympy.S(1)/4)/((c + d*x)**(sympy.S(1)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1718():
    f = 1/((a + b*x)**(sympy.S(7)/4)*(c + d*x)**(sympy.S(5)/4))
    F = -16*d*(a + b*x)**(sympy.S(1)/4)/(3*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2) - 4/((a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(1)/4)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1719():
    f = 1/((a + b*x)**(sympy.S(11)/4)*(c + d*x)**(sympy.S(5)/4))
    F = 128*d**2*(a + b*x)**(sympy.S(1)/4)/(21*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**3) + 32*d/(21*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2) - 4/((a + b*x)**(sympy.S(7)/4)*(c + d*x)**(sympy.S(1)/4)*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1720():
    f = 1/((a + b*x)**(sympy.S(15)/4)*(c + d*x)**(sympy.S(5)/4))
    F = -512*d**3*(a + b*x)**(sympy.S(1)/4)/(77*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**4) - 128*d**2/(77*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**3) + 48*d/(77*(a + b*x)**(sympy.S(7)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2) - 4/((a + b*x)**(sympy.S(11)/4)*(c + d*x)**(sympy.S(1)/4)*(-11*a*d + 11*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1721():
    f = (a + b*x)**(sympy.S(11)/4)/(c + d*x)**(sympy.S(5)/4)
    F = -77*sqrt(2)*b**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(7)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(20*d**(sympy.S(15)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 77*sqrt(2)*b**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(7)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(40*d**(sympy.S(15)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 77*sqrt(b)*sqrt((a + b*x)*(c + d*x))*(-a*d + b*c)*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(10*d**(sympy.S(7)/2)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*(a*d + b*c + 2*b*d*x)) + 22*b*(a + b*x)**(sympy.S(7)/4)*(c + d*x)**(sympy.S(3)/4)/(5*d**2) - 77*b*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)*(-a*d + b*c)/(15*d**3) - 4*(a + b*x)**(sympy.S(11)/4)/(d*(c + d*x)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1722():
    f = (a + b*x)**(sympy.S(7)/4)/(c + d*x)**(sympy.S(5)/4)
    F = 7*sqrt(2)*b**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(5)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(2*d**(sympy.S(11)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 7*sqrt(2)*b**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(5)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(4*d**(sympy.S(11)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 7*sqrt(b)*sqrt((a + b*x)*(c + d*x))*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(d**(sympy.S(5)/2)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*(a*d + b*c + 2*b*d*x)) + 14*b*(a + b*x)**(sympy.S(3)/4)*(c + d*x)**(sympy.S(3)/4)/(3*d**2) - 4*(a + b*x)**(sympy.S(7)/4)/(d*(c + d*x)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1723():
    f = (a + b*x)**(sympy.S(3)/4)/(c + d*x)**(sympy.S(5)/4)
    F = -3*sqrt(2)*b**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(3)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(d**(sympy.S(7)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 3*sqrt(2)*b**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(-a*d + b*c)**(sympy.S(3)/2)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(2*d**(sympy.S(7)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 6*sqrt(b)*sqrt((a + b*x)*(c + d*x))*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(d**(sympy.S(3)/2)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*(a*d + b*c + 2*b*d*x)) - 4*(a + b*x)**(sympy.S(3)/4)/(d*(c + d*x)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1724():
    f = 1/((a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(5)/4))
    F = 2*sqrt(2)*b**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*sqrt(-a*d + b*c)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(d**(sympy.S(3)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - sqrt(2)*b**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*sqrt(-a*d + b*c)*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(d**(sympy.S(3)/4)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 4*sqrt(b)*sqrt((a + b*x)*(c + d*x))*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(sqrt(d)*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*(a*d + b*c + 2*b*d*x)) + 4*(a + b*x)**(sympy.S(3)/4)/((c + d*x)**(sympy.S(1)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1725():
    f = 1/((a + b*x)**(sympy.S(5)/4)*(c + d*x)**(sympy.S(5)/4))
    F = -4*sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/((a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*sqrt(-a*d + b*c)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 2*sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/((a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*sqrt(-a*d + b*c)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) + 8*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/((a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**3*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*(a*d + b*c + 2*b*d*x)) - 8*d*(a + b*x)**(sympy.S(3)/4)/((c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2) - 4/((a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1726():
    f = 1/((a + b*x)**(sympy.S(9)/4)*(c + d*x)**(sympy.S(5)/4))
    F = 24*sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(5)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_e(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(5*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**(sympy.S(3)/2)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 12*sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(5)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)*sqrt((a*d + b*(c + 2*d*x))**2/((-a*d + b*c)**2*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)**2))*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*sqrt((a*d + b*c + 2*b*d*x)**2)*elliptic_f(2*atan(sqrt(2)*b**(sympy.S(1)/4)*d**(sympy.S(1)/4)*((a + b*x)*(c + d*x))**(sympy.S(1)/4)/sqrt(-a*d + b*c)), sympy.S.Half)/(5*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**(sympy.S(3)/2)*(a*d + b*c + 2*b*d*x)*sqrt((a*d + b*(c + 2*d*x))**2)) - 48*sqrt(b)*d**(sympy.S(3)/2)*sqrt((a + b*x)*(c + d*x))*sqrt((a*d + b*(c + 2*d*x))**2)*sqrt((a*d + b*c + 2*b*d*x)**2)/(5*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**4*(2*sqrt(b)*sqrt(d)*sqrt((a + b*x)*(c + d*x))/(-a*d + b*c) + 1)*(a*d + b*c + 2*b*d*x)) + 48*d**2*(a + b*x)**(sympy.S(3)/4)/(5*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**3) + 24*d/(5*(a + b*x)**(sympy.S(1)/4)*(c + d*x)**(sympy.S(1)/4)*(-a*d + b*c)**2) - 4/((a + b*x)**(sympy.S(5)/4)*(c + d*x)**(sympy.S(1)/4)*(-5*a*d + 5*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1727():
    f = 1/((-a*x + 1)**(sympy.S(1)/4)*(b*x + 1)**(sympy.S(3)/4))
    F = -sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(-a*x + 1)**(sympy.S(1)/4)/(b*x + 1)**(sympy.S(1)/4) + sqrt(a) + sqrt(b)*sqrt(-a*x + 1)/sqrt(b*x + 1))/(2*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)) + sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(-a*x + 1)**(sympy.S(1)/4)/(b*x + 1)**(sympy.S(1)/4) + sqrt(a) + sqrt(b)*sqrt(-a*x + 1)/sqrt(b*x + 1))/(2*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)) + sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*(-a*x + 1)**(sympy.S(1)/4)/(a**(sympy.S(1)/4)*(b*x + 1)**(sympy.S(1)/4)))/(a**(sympy.S(1)/4)*b**(sympy.S(3)/4)) - sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*(-a*x + 1)**(sympy.S(1)/4)/(a**(sympy.S(1)/4)*(b*x + 1)**(sympy.S(1)/4)))/(a**(sympy.S(1)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1728():
    f = 1/((-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4))
    F = -sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(2*a) + sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(2*a) - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/a - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1729():
    f = (a + b*x)**(sympy.S(3)/2)/(c + d*x)**(sympy.S(1)/5)
    F = 2*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/5)*(a + b*x)**(sympy.S(5)/2)*hyper((sympy.S(1)/5, sympy.S(5)/2), (sympy.S(7)/2,), -d*(a + b*x)/(-a*d + b*c))/(5*b*(c + d*x)**(sympy.S(1)/5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1730():
    f = sqrt(a + b*x)/(c + d*x)**(sympy.S(1)/5)
    F = 2*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/5)*(a + b*x)**(sympy.S(3)/2)*hyper((sympy.S(1)/5, sympy.S(3)/2), (sympy.S(5)/2,), -d*(a + b*x)/(-a*d + b*c))/(3*b*(c + d*x)**(sympy.S(1)/5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1731():
    f = 1/(sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/5))
    F = 2*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/5)*sqrt(a + b*x)*hyper((sympy.S(1)/5, sympy.S.Half), (sympy.S(3)/2,), -d*(a + b*x)/(-a*d + b*c))/(b*(c + d*x)**(sympy.S(1)/5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1732():
    f = 1/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/5))
    F = -2*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/5)*hyper((sympy.S(-1)/2, sympy.S(1)/5), (sympy.S.Half,), -d*(a + b*x)/(-a*d + b*c))/(b*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1733():
    f = 1/((a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(1)/5))
    F = -2*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/5)*hyper((sympy.S(-3)/2, sympy.S(1)/5), (sympy.S(-1)/2,), -d*(a + b*x)/(-a*d + b*c))/(3*b*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1734():
    f = (a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(1)/6)
    F = 3*(a + b*x)**(sympy.S(7)/2)*(c + d*x)**(sympy.S(1)/6)/(11*b) + (a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(1)/6)*(-3*a*d + 3*b*c)/(176*b*d) - 9*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**2/(352*b*d**2) + 81*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**3/(1408*b*d**3) - 81*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(11)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(2816*b*d**4*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1735():
    f = (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/6)
    F = 3*(a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(1)/6)/(8*b) + (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/6)*(-3*a*d + 3*b*c)/(80*b*d) - 27*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**2/(320*b*d**2) + 27*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(8)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(640*b*d**3*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1736():
    f = sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)
    F = 3*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/6)/(5*b) + sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)*(-3*a*d + 3*b*c)/(20*b*d) - 3*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(5)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(40*b*d**2*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1737():
    f = (c + d*x)**(sympy.S(1)/6)/sqrt(a + b*x)
    F = 3*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)/(2*b) + 3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(4*b*d*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1738():
    f = (c + d*x)**(sympy.S(1)/6)/(a + b*x)**(sympy.S(3)/2)
    F = 3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(3*b*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(1)/3)) - 2*(c + d*x)**(sympy.S(1)/6)/(b*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1739():
    f = (c + d*x)**(sympy.S(1)/6)/(a + b*x)**(sympy.S(5)/2)
    F = -2*3**(sympy.S(3)/4)*d*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(27*b*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(4)/3)) - 2*d*(c + d*x)**(sympy.S(1)/6)/(9*b*sqrt(a + b*x)*(-a*d + b*c)) - 2*(c + d*x)**(sympy.S(1)/6)/(3*b*(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1740():
    f = (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(5)/6)
    F = 3*(a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(5)/6)/(10*b) + (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(5)/6)*(-3*a*d + 3*b*c)/(28*b*d) - 27*sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/6)*(-a*d + b*c)**2/(224*b*d**2) - (81 + 81*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**3/(448*b**(sympy.S(5)/3)*d**2*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) - 81*3**(sympy.S(1)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(10)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(448*b**(sympy.S(5)/3)*d**3*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) - 27*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(10)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(896*b**(sympy.S(5)/3)*d**3*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1741():
    f = sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/6)
    F = 3*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(5)/6)/(7*b) + sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/6)*(-15*a*d + 15*b*c)/(56*b*d) + (45 + 45*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**2/(112*b**(sympy.S(5)/3)*d*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) + 45*3**(sympy.S(1)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(7)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(112*b**(sympy.S(5)/3)*d**2*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) + 15*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(7)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(224*b**(sympy.S(5)/3)*d**2*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1742():
    f = (c + d*x)**(sympy.S(5)/6)/sqrt(a + b*x)
    F = 3*sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/6)/(4*b) - (15 + 15*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)/(8*b**(sympy.S(5)/3)*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) - 15*3**(sympy.S(1)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(4)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(8*b**(sympy.S(5)/3)*d*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) - 5*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(4)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(16*b**(sympy.S(5)/3)*d*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1743():
    f = (c + d*x)**(sympy.S(5)/6)/(a + b*x)**(sympy.S(3)/2)
    F = -2*(c + d*x)**(sympy.S(5)/6)/(b*sqrt(a + b*x)) - d*(5 + 5*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)/(b**(sympy.S(5)/3)*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) - 5*3**(sympy.S(1)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(b**(sympy.S(5)/3)*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) - 3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(5 - 5*sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(6*b**(sympy.S(5)/3)*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1744():
    f = (c + d*x)**(sympy.S(5)/6)/(a + b*x)**(sympy.S(5)/2)
    F = -10*d*(c + d*x)**(sympy.S(5)/6)/(9*b*sqrt(a + b*x)*(-a*d + b*c)) - 2*(c + d*x)**(sympy.S(5)/6)/(3*b*(a + b*x)**(sympy.S(3)/2)) - d**2*(10 + 10*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)/(9*b**(sympy.S(5)/3)*(-a*d + b*c)*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) - 10*3**(sympy.S(1)/4)*d*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(9*b**(sympy.S(5)/3)*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(2)/3)) - 3**(sympy.S(3)/4)*d*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(5 - 5*sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(27*b**(sympy.S(5)/3)*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1745():
    f = (c + d*x)**(sympy.S(5)/6)/(a + b*x)**(sympy.S(7)/2)
    F = 8*d**2*(c + d*x)**(sympy.S(5)/6)/(27*b*sqrt(a + b*x)*(-a*d + b*c)**2) - 2*d*(c + d*x)**(sympy.S(5)/6)/(9*b*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)) - 2*(c + d*x)**(sympy.S(5)/6)/(5*b*(a + b*x)**(sympy.S(5)/2)) + d**3*(8 + 8*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)/(27*b**(sympy.S(5)/3)*(-a*d + b*c)**2*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) + 8*3**(sympy.S(1)/4)*d**2*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(27*b**(sympy.S(5)/3)*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/3)) + 3**(sympy.S(3)/4)*d**2*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(4 - 4*sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(81*b**(sympy.S(5)/3)*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1746():
    f = (a + b*x)**(sympy.S(5)/2)/(c + d*x)**(sympy.S(1)/6)
    F = 3*(a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(5)/6)/(10*d) - (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(5)/6)*(-9*a*d + 9*b*c)/(28*d**2) + 81*sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/6)*(-a*d + b*c)**2/(224*d**3) + (243 + 243*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**3/(448*b**(sympy.S(2)/3)*d**3*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) + 243*3**(sympy.S(1)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(10)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(448*b**(sympy.S(2)/3)*d**4*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) + 81*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(10)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(896*b**(sympy.S(2)/3)*d**4*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1747():
    f = (a + b*x)**(sympy.S(3)/2)/(c + d*x)**(sympy.S(1)/6)
    F = 3*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(5)/6)/(7*d) - sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/6)*(-27*a*d + 27*b*c)/(56*d**2) - (81 + 81*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**2/(112*b**(sympy.S(2)/3)*d**2*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) - 81*3**(sympy.S(1)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(7)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(112*b**(sympy.S(2)/3)*d**3*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) - 27*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(7)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(224*b**(sympy.S(2)/3)*d**3*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1748():
    f = sqrt(a + b*x)/(c + d*x)**(sympy.S(1)/6)
    F = 3*sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/6)/(4*d) + (9 + 9*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)/(8*b**(sympy.S(2)/3)*d*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) + 9*3**(sympy.S(1)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(4)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(8*b**(sympy.S(2)/3)*d**2*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) + 3*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(4)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(16*b**(sympy.S(2)/3)*d**2*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1749():
    f = 1/(sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6))
    F = -(3 + 3*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)/(b**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) - 3*3**(sympy.S(1)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(b**(sympy.S(2)/3)*d*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) - 3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(2*b**(sympy.S(2)/3)*d*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1750():
    f = 1/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/6))
    F = -2*(c + d*x)**(sympy.S(5)/6)/(sqrt(a + b*x)*(-a*d + b*c)) - d*(2 + 2*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)/(b**(sympy.S(2)/3)*(-a*d + b*c)*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) - 2*3**(sympy.S(1)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(b**(sympy.S(2)/3)*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(2)/3)) - 3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(3*b**(sympy.S(2)/3)*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1751():
    f = 1/((a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(1)/6))
    F = 8*d*(c + d*x)**(sympy.S(5)/6)/(9*sqrt(a + b*x)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(5)/6)/((a + b*x)**(sympy.S(3)/2)*(-3*a*d + 3*b*c)) + d**2*(8 + 8*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)/(9*b**(sympy.S(2)/3)*(-a*d + b*c)**2*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) + 8*3**(sympy.S(1)/4)*d*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(9*b**(sympy.S(2)/3)*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/3)) + 3**(sympy.S(3)/4)*d*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(4 - 4*sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(27*b**(sympy.S(2)/3)*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1752():
    f = (a + b*x)**(sympy.S(5)/2)/(c + d*x)**(sympy.S(5)/6)
    F = 3*(a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(1)/6)/(8*d) - (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/6)*(-9*a*d + 9*b*c)/(16*d**2) + 81*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**2/(64*d**3) - 81*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(8)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(128*d**4*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1753():
    f = (a + b*x)**(sympy.S(3)/2)/(c + d*x)**(sympy.S(5)/6)
    F = 3*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/6)/(5*d) - sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)*(-27*a*d + 27*b*c)/(20*d**2) + 27*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(5)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(40*d**3*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1754():
    f = sqrt(a + b*x)/(c + d*x)**(sympy.S(5)/6)
    F = 3*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)/(2*d) - 3*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(4*d**2*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1755():
    f = 1/(sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/6))
    F = 3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(d*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1756():
    f = 1/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(5)/6))
    F = -2*3**(sympy.S(3)/4)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(3*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(4)/3)) - 2*(c + d*x)**(sympy.S(1)/6)/(sqrt(a + b*x)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1757():
    f = 1/((a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(5)/6))
    F = 16*3**(sympy.S(3)/4)*d*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(27*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(7)/3)) + 16*d*(c + d*x)**(sympy.S(1)/6)/(9*sqrt(a + b*x)*(-a*d + b*c)**2) - 2*(c + d*x)**(sympy.S(1)/6)/((a + b*x)**(sympy.S(3)/2)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1758():
    f = (a + b*x)**(sympy.S(5)/2)/(c + d*x)**(sympy.S(7)/6)
    F = -b**(sympy.S(1)/3)*(1215 + 1215*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**2/(112*d**3*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) - 1215*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(7)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(112*d**4*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) - 405*3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(7)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(224*d**4*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) + 45*b*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(5)/6)/(7*d**2) - 405*b*sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/6)*(-a*d + b*c)/(56*d**3) - 6*(a + b*x)**(sympy.S(5)/2)/(d*(c + d*x)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1759():
    f = (a + b*x)**(sympy.S(3)/2)/(c + d*x)**(sympy.S(7)/6)
    F = b**(sympy.S(1)/3)*(81 + 81*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)/(8*d**2*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) + 81*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(4)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(8*d**3*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) + 27*3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(4)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(16*d**3*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) + 27*b*sqrt(a + b*x)*(c + d*x)**(sympy.S(5)/6)/(4*d**2) - 6*(a + b*x)**(sympy.S(3)/2)/(d*(c + d*x)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1760():
    f = sqrt(a + b*x)/(c + d*x)**(sympy.S(7)/6)
    F = -b**(sympy.S(1)/3)*(9 + 9*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)/(d*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) - 9*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(d**2*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) - 3*3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(2*d**2*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)) - 6*sqrt(a + b*x)/(d*(c + d*x)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1761():
    f = 1/(sqrt(a + b*x)*(c + d*x)**(sympy.S(7)/6))
    F = b**(sympy.S(1)/3)*(6 + 6*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)/((-a*d + b*c)*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) + 6*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(d*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(2)/3)) + 3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(d*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(2)/3)) + 6*sqrt(a + b*x)/((c + d*x)**(sympy.S(1)/6)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1762():
    f = 1/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(7)/6))
    F = -b**(sympy.S(1)/3)*d*(8 + 8*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)/((-a*d + b*c)**2*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) - 8*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/3)) - 3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(4 - 4*sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(3*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(5)/3)) - 8*d*sqrt(a + b*x)/((c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**2) - 2/(sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1763():
    f = 1/((a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(7)/6))
    F = b**(sympy.S(1)/3)*d**2*(80 + 80*sqrt(3))*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)/(9*(-a*d + b*c)**3*(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))) + 80*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*d*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_e(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(9*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(8)/3)) + 3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*d*sqrt((b**(sympy.S(2)/3)*(c + d*x)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-a*d + b*c)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*(40 - 40*sqrt(3))*(c + d*x)**(sympy.S(1)/6)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(1 - sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(27*sqrt(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))/(-b**(sympy.S(1)/3)*(1 + sqrt(3))*(c + d*x)**(sympy.S(1)/3) + (-a*d + b*c)**(sympy.S(1)/3))**2)*sqrt(a + b*x)*(-a*d + b*c)**(sympy.S(8)/3)) + 80*d**2*sqrt(a + b*x)/(9*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**3) + 20*d/(9*sqrt(a + b*x)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**2) - 2/((a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(1)/6)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1764():
    f = (a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(13)/6)
    F = 6*(a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**2*hyper((sympy.S(-13)/6, sympy.S(7)/6), (sympy.S(13)/6,), -d*(a + b*x)/(-a*d + b*c))/(7*b**3*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1765():
    f = (a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(7)/6)
    F = (a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(1)/6)*(-6*a*d + 6*b*c)*hyper((sympy.S(-7)/6, sympy.S(7)/6), (sympy.S(13)/6,), -d*(a + b*x)/(-a*d + b*c))/(7*b**2*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1766():
    f = (a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)
    F = 6*(a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(1)/6)*hyper((sympy.S(-1)/6, sympy.S(7)/6), (sympy.S(13)/6,), -d*(a + b*x)/(-a*d + b*c))/(7*b*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1767():
    f = (a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(5)/6)
    F = 6*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(5)/6)*(a + b*x)**(sympy.S(7)/6)*hyper((sympy.S(5)/6, sympy.S(7)/6), (sympy.S(13)/6,), -d*(a + b*x)/(-a*d + b*c))/(7*b*(c + d*x)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1768():
    f = (a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(11)/6)
    F = 6*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(5)/6)*(a + b*x)**(sympy.S(7)/6)*hyper((sympy.S(7)/6, sympy.S(11)/6), (sympy.S(13)/6,), -d*(a + b*x)/(-a*d + b*c))/((c + d*x)**(sympy.S(5)/6)*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1769():
    f = (a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(17)/6)
    F = 6*b*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(5)/6)*(a + b*x)**(sympy.S(7)/6)*hyper((sympy.S(7)/6, sympy.S(17)/6), (sympy.S(13)/6,), -d*(a + b*x)/(-a*d + b*c))/(7*(c + d*x)**(sympy.S(5)/6)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1770():
    f = (a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(5)/6)
    F = (a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(5)/6)/(2*b) + (a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(5)/6)*(-5*a*d + 5*b*c)/(12*b*d) + 5*(-a*d + b*c)**2*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(144*b**(sympy.S(11)/6)*d**(sympy.S(7)/6)) - 5*(-a*d + b*c)**2*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(144*b**(sympy.S(11)/6)*d**(sympy.S(7)/6)) + 5*sqrt(3)*(-a*d + b*c)**2*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(72*b**(sympy.S(11)/6)*d**(sympy.S(7)/6)) - 5*sqrt(3)*(-a*d + b*c)**2*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(72*b**(sympy.S(11)/6)*d**(sympy.S(7)/6)) - 5*(-a*d + b*c)**2*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(36*b**(sympy.S(11)/6)*d**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1771():
    f = (a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6)
    F = (a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(5)/6)/d + (-a*d + b*c)*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(12*b**(sympy.S(5)/6)*d**(sympy.S(7)/6)) - (-a*d + b*c)*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(12*b**(sympy.S(5)/6)*d**(sympy.S(7)/6)) + sqrt(3)*(-a*d + b*c)*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(6*b**(sympy.S(5)/6)*d**(sympy.S(7)/6)) - sqrt(3)*(-a*d + b*c)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(6*b**(sympy.S(5)/6)*d**(sympy.S(7)/6)) - (-a*d + b*c)*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(3*b**(sympy.S(5)/6)*d**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1772():
    f = (a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(7)/6)
    F = -b**(sympy.S(1)/6)*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(2*d**(sympy.S(7)/6)) + b**(sympy.S(1)/6)*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(2*d**(sympy.S(7)/6)) - sqrt(3)*b**(sympy.S(1)/6)*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/d**(sympy.S(7)/6) + sqrt(3)*b**(sympy.S(1)/6)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/d**(sympy.S(7)/6) + 2*b**(sympy.S(1)/6)*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/d**(sympy.S(7)/6) - 6*(a + b*x)**(sympy.S(1)/6)/(d*(c + d*x)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1773():
    f = (a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(13)/6)
    F = 6*(a + b*x)**(sympy.S(7)/6)/((c + d*x)**(sympy.S(7)/6)*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1774():
    f = (a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(19)/6)
    F = 36*b*(a + b*x)**(sympy.S(7)/6)/(91*(c + d*x)**(sympy.S(7)/6)*(-a*d + b*c)**2) + 6*(a + b*x)**(sympy.S(7)/6)/((c + d*x)**(sympy.S(13)/6)*(-13*a*d + 13*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1775():
    f = (a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(25)/6)
    F = 432*b**2*(a + b*x)**(sympy.S(7)/6)/(1729*(c + d*x)**(sympy.S(7)/6)*(-a*d + b*c)**3) + 72*b*(a + b*x)**(sympy.S(7)/6)/(247*(c + d*x)**(sympy.S(13)/6)*(-a*d + b*c)**2) + 6*(a + b*x)**(sympy.S(7)/6)/((c + d*x)**(sympy.S(19)/6)*(-19*a*d + 19*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1776():
    f = (a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(31)/6)
    F = 7776*b**3*(a + b*x)**(sympy.S(7)/6)/(43225*(c + d*x)**(sympy.S(7)/6)*(-a*d + b*c)**4) + 1296*b**2*(a + b*x)**(sympy.S(7)/6)/(6175*(c + d*x)**(sympy.S(13)/6)*(-a*d + b*c)**3) + 108*b*(a + b*x)**(sympy.S(7)/6)/(475*(c + d*x)**(sympy.S(19)/6)*(-a*d + b*c)**2) + 6*(a + b*x)**(sympy.S(7)/6)/((c + d*x)**(sympy.S(25)/6)*(-25*a*d + 25*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1777():
    f = (a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(1)/6)
    F = (a + b*x)**(sympy.S(11)/6)*(c + d*x)**(sympy.S(1)/6)/(2*b) + (a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)/(12*b*d) + 5*(-a*d + b*c)**2*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(144*b**(sympy.S(7)/6)*d**(sympy.S(11)/6)) - 5*(-a*d + b*c)**2*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(144*b**(sympy.S(7)/6)*d**(sympy.S(11)/6)) - 5*sqrt(3)*(-a*d + b*c)**2*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(72*b**(sympy.S(7)/6)*d**(sympy.S(11)/6)) + 5*sqrt(3)*(-a*d + b*c)**2*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(72*b**(sympy.S(7)/6)*d**(sympy.S(11)/6)) - 5*(-a*d + b*c)**2*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(36*b**(sympy.S(7)/6)*d**(sympy.S(11)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1778():
    f = (a + b*x)**(sympy.S(5)/6)/(c + d*x)**(sympy.S(5)/6)
    F = (a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(1)/6)/d + (-5*a*d + 5*b*c)*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(12*b**(sympy.S(1)/6)*d**(sympy.S(11)/6)) - (-5*a*d + 5*b*c)*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(12*b**(sympy.S(1)/6)*d**(sympy.S(11)/6)) - sqrt(3)*(-5*a*d + 5*b*c)*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(6*b**(sympy.S(1)/6)*d**(sympy.S(11)/6)) + sqrt(3)*(-5*a*d + 5*b*c)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(6*b**(sympy.S(1)/6)*d**(sympy.S(11)/6)) - (-5*a*d + 5*b*c)*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(3*b**(sympy.S(1)/6)*d**(sympy.S(11)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1779():
    f = (a + b*x)**(sympy.S(5)/6)/(c + d*x)**(sympy.S(11)/6)
    F = -b**(sympy.S(5)/6)*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(2*d**(sympy.S(11)/6)) + b**(sympy.S(5)/6)*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(2*d**(sympy.S(11)/6)) + sqrt(3)*b**(sympy.S(5)/6)*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/d**(sympy.S(11)/6) - sqrt(3)*b**(sympy.S(5)/6)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/d**(sympy.S(11)/6) + 2*b**(sympy.S(5)/6)*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/d**(sympy.S(11)/6) - 6*(a + b*x)**(sympy.S(5)/6)/(5*d*(c + d*x)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1780():
    f = (a + b*x)**(sympy.S(5)/6)/(c + d*x)**(sympy.S(17)/6)
    F = 6*(a + b*x)**(sympy.S(11)/6)/((c + d*x)**(sympy.S(11)/6)*(-11*a*d + 11*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1781():
    f = (a + b*x)**(sympy.S(5)/6)/(c + d*x)**(sympy.S(23)/6)
    F = 36*b*(a + b*x)**(sympy.S(11)/6)/(187*(c + d*x)**(sympy.S(11)/6)*(-a*d + b*c)**2) + 6*(a + b*x)**(sympy.S(11)/6)/((c + d*x)**(sympy.S(17)/6)*(-17*a*d + 17*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1782():
    f = (a + b*x)**(sympy.S(5)/6)/(c + d*x)**(sympy.S(29)/6)
    F = 432*b**2*(a + b*x)**(sympy.S(11)/6)/(4301*(c + d*x)**(sympy.S(11)/6)*(-a*d + b*c)**3) + 72*b*(a + b*x)**(sympy.S(11)/6)/(391*(c + d*x)**(sympy.S(17)/6)*(-a*d + b*c)**2) + 6*(a + b*x)**(sympy.S(11)/6)/((c + d*x)**(sympy.S(23)/6)*(-23*a*d + 23*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1783():
    f = (a + b*x)**(sympy.S(5)/6)/(c + d*x)**(sympy.S(35)/6)
    F = 7776*b**3*(a + b*x)**(sympy.S(11)/6)/(124729*(c + d*x)**(sympy.S(11)/6)*(-a*d + b*c)**4) + 1296*b**2*(a + b*x)**(sympy.S(11)/6)/(11339*(c + d*x)**(sympy.S(17)/6)*(-a*d + b*c)**3) + 108*b*(a + b*x)**(sympy.S(11)/6)/(667*(c + d*x)**(sympy.S(23)/6)*(-a*d + b*c)**2) + 6*(a + b*x)**(sympy.S(11)/6)/((c + d*x)**(sympy.S(29)/6)*(-29*a*d + 29*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1784():
    f = (a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(11)/6)
    F = (a + b*x)**(sympy.S(11)/6)*(c + d*x)**(sympy.S(5)/6)*(-6*a*d + 6*b*c)*hyper((sympy.S(-11)/6, sympy.S(11)/6), (sympy.S(17)/6,), -d*(a + b*x)/(-a*d + b*c))/(11*b**2*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1785():
    f = (a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(5)/6)
    F = 6*(a + b*x)**(sympy.S(11)/6)*(c + d*x)**(sympy.S(5)/6)*hyper((sympy.S(-5)/6, sympy.S(11)/6), (sympy.S(17)/6,), -d*(a + b*x)/(-a*d + b*c))/(11*b*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1786():
    f = (a + b*x)**(sympy.S(5)/6)/(c + d*x)**(sympy.S(1)/6)
    F = 6*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6)*(a + b*x)**(sympy.S(11)/6)*hyper((sympy.S(1)/6, sympy.S(11)/6), (sympy.S(17)/6,), -d*(a + b*x)/(-a*d + b*c))/(11*b*(c + d*x)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1787():
    f = (a + b*x)**(sympy.S(5)/6)/(c + d*x)**(sympy.S(7)/6)
    F = 6*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6)*(a + b*x)**(sympy.S(11)/6)*hyper((sympy.S(7)/6, sympy.S(11)/6), (sympy.S(17)/6,), -d*(a + b*x)/(-a*d + b*c))/((c + d*x)**(sympy.S(1)/6)*(-11*a*d + 11*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1788():
    f = (a + b*x)**(sympy.S(5)/6)/(c + d*x)**(sympy.S(13)/6)
    F = 6*b*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6)*(a + b*x)**(sympy.S(11)/6)*hyper((sympy.S(11)/6, sympy.S(13)/6), (sympy.S(17)/6,), -d*(a + b*x)/(-a*d + b*c))/(11*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1789():
    f = (a + b*x)**(sympy.S(5)/6)/(c + d*x)**(sympy.S(19)/6)
    F = 6*b**2*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6)*(a + b*x)**(sympy.S(11)/6)*hyper((sympy.S(11)/6, sympy.S(19)/6), (sympy.S(17)/6,), -d*(a + b*x)/(-a*d + b*c))/(11*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1790():
    f = (a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(13)/6)
    F = 6*(a + b*x)**(sympy.S(13)/6)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**2*hyper((sympy.S(-13)/6, sympy.S(13)/6), (sympy.S(19)/6,), -d*(a + b*x)/(-a*d + b*c))/(13*b**3*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1791():
    f = (a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(7)/6)
    F = (a + b*x)**(sympy.S(13)/6)*(c + d*x)**(sympy.S(1)/6)*(-6*a*d + 6*b*c)*hyper((sympy.S(-7)/6, sympy.S(13)/6), (sympy.S(19)/6,), -d*(a + b*x)/(-a*d + b*c))/(13*b**2*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1792():
    f = (a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(1)/6)
    F = 6*(a + b*x)**(sympy.S(13)/6)*(c + d*x)**(sympy.S(1)/6)*hyper((sympy.S(-1)/6, sympy.S(13)/6), (sympy.S(19)/6,), -d*(a + b*x)/(-a*d + b*c))/(13*b*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1793():
    f = (a + b*x)**(sympy.S(7)/6)/(c + d*x)**(sympy.S(5)/6)
    F = 6*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(5)/6)*(a + b*x)**(sympy.S(13)/6)*hyper((sympy.S(5)/6, sympy.S(13)/6), (sympy.S(19)/6,), -d*(a + b*x)/(-a*d + b*c))/(13*b*(c + d*x)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1794():
    f = (a + b*x)**(sympy.S(7)/6)/(c + d*x)**(sympy.S(11)/6)
    F = 6*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(5)/6)*(a + b*x)**(sympy.S(13)/6)*hyper((sympy.S(11)/6, sympy.S(13)/6), (sympy.S(19)/6,), -d*(a + b*x)/(-a*d + b*c))/((c + d*x)**(sympy.S(5)/6)*(-13*a*d + 13*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1795():
    f = (a + b*x)**(sympy.S(7)/6)/(c + d*x)**(sympy.S(17)/6)
    F = 6*b*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(5)/6)*(a + b*x)**(sympy.S(13)/6)*hyper((sympy.S(13)/6, sympy.S(17)/6), (sympy.S(19)/6,), -d*(a + b*x)/(-a*d + b*c))/(13*(c + d*x)**(sympy.S(5)/6)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1796():
    f = (a + b*x)**(sympy.S(7)/6)/(c + d*x)**(sympy.S(1)/6)
    F = (a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(5)/6)/(2*d) - (a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(5)/6)*(-7*a*d + 7*b*c)/(12*d**2) - 7*(-a*d + b*c)**2*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(144*b**(sympy.S(5)/6)*d**(sympy.S(13)/6)) + 7*(-a*d + b*c)**2*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(144*b**(sympy.S(5)/6)*d**(sympy.S(13)/6)) - 7*sqrt(3)*(-a*d + b*c)**2*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(72*b**(sympy.S(5)/6)*d**(sympy.S(13)/6)) + 7*sqrt(3)*(-a*d + b*c)**2*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(72*b**(sympy.S(5)/6)*d**(sympy.S(13)/6)) + 7*(-a*d + b*c)**2*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(36*b**(sympy.S(5)/6)*d**(sympy.S(13)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1797():
    f = (a + b*x)**(sympy.S(7)/6)/(c + d*x)**(sympy.S(7)/6)
    F = 7*b**(sympy.S(1)/6)*(-a*d + b*c)*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(12*d**(sympy.S(13)/6)) - 7*b**(sympy.S(1)/6)*(-a*d + b*c)*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(12*d**(sympy.S(13)/6)) + 7*sqrt(3)*b**(sympy.S(1)/6)*(-a*d + b*c)*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(6*d**(sympy.S(13)/6)) - 7*sqrt(3)*b**(sympy.S(1)/6)*(-a*d + b*c)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(6*d**(sympy.S(13)/6)) - 7*b**(sympy.S(1)/6)*(-a*d + b*c)*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(3*d**(sympy.S(13)/6)) + 7*b*(a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(5)/6)/d**2 - 6*(a + b*x)**(sympy.S(7)/6)/(d*(c + d*x)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1798():
    f = (a + b*x)**(sympy.S(7)/6)/(c + d*x)**(sympy.S(13)/6)
    F = -b**(sympy.S(7)/6)*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(2*d**(sympy.S(13)/6)) + b**(sympy.S(7)/6)*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(2*d**(sympy.S(13)/6)) - sqrt(3)*b**(sympy.S(7)/6)*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/d**(sympy.S(13)/6) + sqrt(3)*b**(sympy.S(7)/6)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/d**(sympy.S(13)/6) + 2*b**(sympy.S(7)/6)*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/d**(sympy.S(13)/6) - 6*b*(a + b*x)**(sympy.S(1)/6)/(d**2*(c + d*x)**(sympy.S(1)/6)) - 6*(a + b*x)**(sympy.S(7)/6)/(7*d*(c + d*x)**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1799():
    f = (a + b*x)**(sympy.S(7)/6)/(c + d*x)**(sympy.S(19)/6)
    F = 6*(a + b*x)**(sympy.S(13)/6)/((c + d*x)**(sympy.S(13)/6)*(-13*a*d + 13*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1800():
    f = (a + b*x)**(sympy.S(7)/6)/(c + d*x)**(sympy.S(25)/6)
    F = 36*b*(a + b*x)**(sympy.S(13)/6)/(247*(c + d*x)**(sympy.S(13)/6)*(-a*d + b*c)**2) + 6*(a + b*x)**(sympy.S(13)/6)/((c + d*x)**(sympy.S(19)/6)*(-19*a*d + 19*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1801():
    f = (a + b*x)**(sympy.S(7)/6)/(c + d*x)**(sympy.S(31)/6)
    F = 432*b**2*(a + b*x)**(sympy.S(13)/6)/(6175*(c + d*x)**(sympy.S(13)/6)*(-a*d + b*c)**3) + 72*b*(a + b*x)**(sympy.S(13)/6)/(475*(c + d*x)**(sympy.S(19)/6)*(-a*d + b*c)**2) + 6*(a + b*x)**(sympy.S(13)/6)/((c + d*x)**(sympy.S(25)/6)*(-25*a*d + 25*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1802():
    f = (a + b*x)**(sympy.S(7)/6)/(c + d*x)**(sympy.S(37)/6)
    F = 7776*b**3*(a + b*x)**(sympy.S(13)/6)/(191425*(c + d*x)**(sympy.S(13)/6)*(-a*d + b*c)**4) + 1296*b**2*(a + b*x)**(sympy.S(13)/6)/(14725*(c + d*x)**(sympy.S(19)/6)*(-a*d + b*c)**3) + 108*b*(a + b*x)**(sympy.S(13)/6)/(775*(c + d*x)**(sympy.S(25)/6)*(-a*d + b*c)**2) + 6*(a + b*x)**(sympy.S(13)/6)/((c + d*x)**(sympy.S(31)/6)*(-31*a*d + 31*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1803():
    f = (c + d*x)**(sympy.S(7)/6)/(a + b*x)**(sympy.S(1)/6)
    F = (a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(7)/6)/(2*b) + (a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(1)/6)*(-7*a*d + 7*b*c)/(12*b**2) - 7*(-a*d + b*c)**2*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(144*b**(sympy.S(13)/6)*d**(sympy.S(5)/6)) + 7*(-a*d + b*c)**2*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(144*b**(sympy.S(13)/6)*d**(sympy.S(5)/6)) + 7*sqrt(3)*(-a*d + b*c)**2*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(72*b**(sympy.S(13)/6)*d**(sympy.S(5)/6)) - 7*sqrt(3)*(-a*d + b*c)**2*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(72*b**(sympy.S(13)/6)*d**(sympy.S(5)/6)) + 7*(-a*d + b*c)**2*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(36*b**(sympy.S(13)/6)*d**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1804():
    f = (c + d*x)**(sympy.S(1)/6)/(a + b*x)**(sympy.S(1)/6)
    F = (a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(1)/6)/b - (-a*d + b*c)*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(12*b**(sympy.S(7)/6)*d**(sympy.S(5)/6)) + (-a*d + b*c)*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(12*b**(sympy.S(7)/6)*d**(sympy.S(5)/6)) + sqrt(3)*(-a*d + b*c)*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(6*b**(sympy.S(7)/6)*d**(sympy.S(5)/6)) - sqrt(3)*(-a*d + b*c)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(6*b**(sympy.S(7)/6)*d**(sympy.S(5)/6)) + (-a*d + b*c)*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(3*b**(sympy.S(7)/6)*d**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1805():
    f = 1/((a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(5)/6))
    F = -log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(2*b**(sympy.S(1)/6)*d**(sympy.S(5)/6)) + log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(2*b**(sympy.S(1)/6)*d**(sympy.S(5)/6)) + sqrt(3)*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(b**(sympy.S(1)/6)*d**(sympy.S(5)/6)) - sqrt(3)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(b**(sympy.S(1)/6)*d**(sympy.S(5)/6)) + 2*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(b**(sympy.S(1)/6)*d**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1806():
    f = 1/((a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(11)/6))
    F = 6*(a + b*x)**(sympy.S(5)/6)/((c + d*x)**(sympy.S(5)/6)*(-5*a*d + 5*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1807():
    f = 1/((a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(17)/6))
    F = 36*b*(a + b*x)**(sympy.S(5)/6)/(55*(c + d*x)**(sympy.S(5)/6)*(-a*d + b*c)**2) + 6*(a + b*x)**(sympy.S(5)/6)/((c + d*x)**(sympy.S(11)/6)*(-11*a*d + 11*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1808():
    f = 1/((a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(23)/6))
    F = 432*b**2*(a + b*x)**(sympy.S(5)/6)/(935*(c + d*x)**(sympy.S(5)/6)*(-a*d + b*c)**3) + 72*b*(a + b*x)**(sympy.S(5)/6)/(187*(c + d*x)**(sympy.S(11)/6)*(-a*d + b*c)**2) + 6*(a + b*x)**(sympy.S(5)/6)/((c + d*x)**(sympy.S(17)/6)*(-17*a*d + 17*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1809():
    f = 1/((a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(29)/6))
    F = 7776*b**3*(a + b*x)**(sympy.S(5)/6)/(21505*(c + d*x)**(sympy.S(5)/6)*(-a*d + b*c)**4) + 1296*b**2*(a + b*x)**(sympy.S(5)/6)/(4301*(c + d*x)**(sympy.S(11)/6)*(-a*d + b*c)**3) + 108*b*(a + b*x)**(sympy.S(5)/6)/(391*(c + d*x)**(sympy.S(17)/6)*(-a*d + b*c)**2) + 6*(a + b*x)**(sympy.S(5)/6)/((c + d*x)**(sympy.S(23)/6)*(-23*a*d + 23*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1810():
    f = (c + d*x)**(sympy.S(11)/6)/(a + b*x)**(sympy.S(1)/6)
    F = (a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(5)/6)*(-6*a*d + 6*b*c)*hyper((sympy.S(-11)/6, sympy.S(5)/6), (sympy.S(11)/6,), -d*(a + b*x)/(-a*d + b*c))/(5*b**2*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1811():
    f = (c + d*x)**(sympy.S(5)/6)/(a + b*x)**(sympy.S(1)/6)
    F = 6*(a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(5)/6)*hyper((sympy.S(-5)/6, sympy.S(5)/6), (sympy.S(11)/6,), -d*(a + b*x)/(-a*d + b*c))/(5*b*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1812():
    f = 1/((a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6))
    F = 6*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6)*(a + b*x)**(sympy.S(5)/6)*hyper((sympy.S(1)/6, sympy.S(5)/6), (sympy.S(11)/6,), -d*(a + b*x)/(-a*d + b*c))/(5*b*(c + d*x)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1813():
    f = 1/((a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(7)/6))
    F = 6*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6)*(a + b*x)**(sympy.S(5)/6)*hyper((sympy.S(5)/6, sympy.S(7)/6), (sympy.S(11)/6,), -d*(a + b*x)/(-a*d + b*c))/((c + d*x)**(sympy.S(1)/6)*(-5*a*d + 5*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1814():
    f = 1/((a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(13)/6))
    F = 6*b*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6)*(a + b*x)**(sympy.S(5)/6)*hyper((sympy.S(5)/6, sympy.S(13)/6), (sympy.S(11)/6,), -d*(a + b*x)/(-a*d + b*c))/(5*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1815():
    f = 1/((a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(19)/6))
    F = 6*b**2*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6)*(a + b*x)**(sympy.S(5)/6)*hyper((sympy.S(5)/6, sympy.S(19)/6), (sympy.S(11)/6,), -d*(a + b*x)/(-a*d + b*c))/(5*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1816():
    f = (c + d*x)**(sympy.S(13)/6)/(a + b*x)**(sympy.S(5)/6)
    F = 6*(a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**2*hyper((sympy.S(-13)/6, sympy.S(1)/6), (sympy.S(7)/6,), -d*(a + b*x)/(-a*d + b*c))/(b**3*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1817():
    f = (c + d*x)**(sympy.S(7)/6)/(a + b*x)**(sympy.S(5)/6)
    F = (a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)*(-6*a*d + 6*b*c)*hyper((sympy.S(-7)/6, sympy.S(1)/6), (sympy.S(7)/6,), -d*(a + b*x)/(-a*d + b*c))/(b**2*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1818():
    f = (c + d*x)**(sympy.S(1)/6)/(a + b*x)**(sympy.S(5)/6)
    F = 6*(a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)*hyper((sympy.S(-1)/6, sympy.S(1)/6), (sympy.S(7)/6,), -d*(a + b*x)/(-a*d + b*c))/(b*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1819():
    f = 1/((a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(5)/6))
    F = 6*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(5)/6)*(a + b*x)**(sympy.S(1)/6)*hyper((sympy.S(1)/6, sympy.S(5)/6), (sympy.S(7)/6,), -d*(a + b*x)/(-a*d + b*c))/(b*(c + d*x)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1820():
    f = 1/((a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(11)/6))
    F = 6*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(5)/6)*(a + b*x)**(sympy.S(1)/6)*hyper((sympy.S(1)/6, sympy.S(11)/6), (sympy.S(7)/6,), -d*(a + b*x)/(-a*d + b*c))/((c + d*x)**(sympy.S(5)/6)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1821():
    f = 1/((a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(17)/6))
    F = 6*b*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(5)/6)*(a + b*x)**(sympy.S(1)/6)*hyper((sympy.S(1)/6, sympy.S(17)/6), (sympy.S(7)/6,), -d*(a + b*x)/(-a*d + b*c))/((c + d*x)**(sympy.S(5)/6)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1822():
    f = (c + d*x)**(sympy.S(11)/6)/(a + b*x)**(sympy.S(5)/6)
    F = (a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(11)/6)/(2*b) + (a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(5)/6)*(-11*a*d + 11*b*c)/(12*b**2) - 55*(-a*d + b*c)**2*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(144*b**(sympy.S(17)/6)*d**(sympy.S(1)/6)) + 55*(-a*d + b*c)**2*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(144*b**(sympy.S(17)/6)*d**(sympy.S(1)/6)) - 55*sqrt(3)*(-a*d + b*c)**2*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(72*b**(sympy.S(17)/6)*d**(sympy.S(1)/6)) + 55*sqrt(3)*(-a*d + b*c)**2*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(72*b**(sympy.S(17)/6)*d**(sympy.S(1)/6)) + 55*(-a*d + b*c)**2*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(36*b**(sympy.S(17)/6)*d**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1823():
    f = (c + d*x)**(sympy.S(5)/6)/(a + b*x)**(sympy.S(5)/6)
    F = (a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(5)/6)/b - (-5*a*d + 5*b*c)*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(12*b**(sympy.S(11)/6)*d**(sympy.S(1)/6)) + (-5*a*d + 5*b*c)*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(12*b**(sympy.S(11)/6)*d**(sympy.S(1)/6)) - sqrt(3)*(-5*a*d + 5*b*c)*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(6*b**(sympy.S(11)/6)*d**(sympy.S(1)/6)) + sqrt(3)*(-5*a*d + 5*b*c)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(6*b**(sympy.S(11)/6)*d**(sympy.S(1)/6)) + (-5*a*d + 5*b*c)*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(3*b**(sympy.S(11)/6)*d**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1824():
    f = 1/((a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(1)/6))
    F = -log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(2*b**(sympy.S(5)/6)*d**(sympy.S(1)/6)) + log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(2*b**(sympy.S(5)/6)*d**(sympy.S(1)/6)) - sqrt(3)*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(b**(sympy.S(5)/6)*d**(sympy.S(1)/6)) + sqrt(3)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(b**(sympy.S(5)/6)*d**(sympy.S(1)/6)) + 2*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(b**(sympy.S(5)/6)*d**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1825():
    f = 1/((a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(7)/6))
    F = 6*(a + b*x)**(sympy.S(1)/6)/((c + d*x)**(sympy.S(1)/6)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1826():
    f = 1/((a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(13)/6))
    F = 36*b*(a + b*x)**(sympy.S(1)/6)/(7*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**2) + 6*(a + b*x)**(sympy.S(1)/6)/((c + d*x)**(sympy.S(7)/6)*(-7*a*d + 7*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1827():
    f = 1/((a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(19)/6))
    F = 432*b**2*(a + b*x)**(sympy.S(1)/6)/(91*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**3) + 72*b*(a + b*x)**(sympy.S(1)/6)/(91*(c + d*x)**(sympy.S(7)/6)*(-a*d + b*c)**2) + 6*(a + b*x)**(sympy.S(1)/6)/((c + d*x)**(sympy.S(13)/6)*(-13*a*d + 13*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1828():
    f = 1/((a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(25)/6))
    F = 7776*b**3*(a + b*x)**(sympy.S(1)/6)/(1729*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**4) + 1296*b**2*(a + b*x)**(sympy.S(1)/6)/(1729*(c + d*x)**(sympy.S(7)/6)*(-a*d + b*c)**3) + 108*b*(a + b*x)**(sympy.S(1)/6)/(247*(c + d*x)**(sympy.S(13)/6)*(-a*d + b*c)**2) + 6*(a + b*x)**(sympy.S(1)/6)/((c + d*x)**(sympy.S(19)/6)*(-19*a*d + 19*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1829():
    f = (c + d*x)**(sympy.S(13)/6)/(a + b*x)**(sympy.S(7)/6)
    F = -6*(c + d*x)**(sympy.S(13)/6)/(b*(a + b*x)**(sympy.S(1)/6)) + 13*d*(a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(7)/6)/(2*b**2) + 91*d*(a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)/(12*b**3) - 91*d**(sympy.S(1)/6)*(-a*d + b*c)**2*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(144*b**(sympy.S(19)/6)) + 91*d**(sympy.S(1)/6)*(-a*d + b*c)**2*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(144*b**(sympy.S(19)/6)) + 91*sqrt(3)*d**(sympy.S(1)/6)*(-a*d + b*c)**2*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(72*b**(sympy.S(19)/6)) - 91*sqrt(3)*d**(sympy.S(1)/6)*(-a*d + b*c)**2*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(72*b**(sympy.S(19)/6)) + 91*d**(sympy.S(1)/6)*(-a*d + b*c)**2*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(36*b**(sympy.S(19)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1830():
    f = (c + d*x)**(sympy.S(7)/6)/(a + b*x)**(sympy.S(7)/6)
    F = -6*(c + d*x)**(sympy.S(7)/6)/(b*(a + b*x)**(sympy.S(1)/6)) + 7*d*(a + b*x)**(sympy.S(5)/6)*(c + d*x)**(sympy.S(1)/6)/b**2 - 7*d**(sympy.S(1)/6)*(-a*d + b*c)*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(12*b**(sympy.S(13)/6)) + 7*d**(sympy.S(1)/6)*(-a*d + b*c)*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(12*b**(sympy.S(13)/6)) + 7*sqrt(3)*d**(sympy.S(1)/6)*(-a*d + b*c)*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(6*b**(sympy.S(13)/6)) - 7*sqrt(3)*d**(sympy.S(1)/6)*(-a*d + b*c)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(6*b**(sympy.S(13)/6)) + 7*d**(sympy.S(1)/6)*(-a*d + b*c)*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/(3*b**(sympy.S(13)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1831():
    f = (c + d*x)**(sympy.S(1)/6)/(a + b*x)**(sympy.S(7)/6)
    F = -6*(c + d*x)**(sympy.S(1)/6)/(b*(a + b*x)**(sympy.S(1)/6)) - d**(sympy.S(1)/6)*log(-b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(2*b**(sympy.S(7)/6)) + d**(sympy.S(1)/6)*log(b**(sympy.S(1)/6)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(c + d*x)**(sympy.S(1)/6) + b**(sympy.S(1)/3) + d**(sympy.S(1)/3)*(a + b*x)**(sympy.S(1)/3)/(c + d*x)**(sympy.S(1)/3))/(2*b**(sympy.S(7)/6)) + sqrt(3)*d**(sympy.S(1)/6)*atan(sqrt(3)/3 - 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/b**(sympy.S(7)/6) - sqrt(3)*d**(sympy.S(1)/6)*atan(sqrt(3)/3 + 2*sqrt(3)*d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(3*b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/b**(sympy.S(7)/6) + 2*d**(sympy.S(1)/6)*atanh(d**(sympy.S(1)/6)*(a + b*x)**(sympy.S(1)/6)/(b**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)))/b**(sympy.S(7)/6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1832():
    f = 1/((a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(5)/6))
    F = -6*(c + d*x)**(sympy.S(1)/6)/((a + b*x)**(sympy.S(1)/6)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1833():
    f = 1/((a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(11)/6))
    F = -36*d*(a + b*x)**(sympy.S(5)/6)/(5*(c + d*x)**(sympy.S(5)/6)*(-a*d + b*c)**2) - 6/((a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(5)/6)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1834():
    f = 1/((a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(17)/6))
    F = -432*b*d*(a + b*x)**(sympy.S(5)/6)/(55*(c + d*x)**(sympy.S(5)/6)*(-a*d + b*c)**3) - 72*d*(a + b*x)**(sympy.S(5)/6)/(11*(c + d*x)**(sympy.S(11)/6)*(-a*d + b*c)**2) - 6/((a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(11)/6)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1835():
    f = 1/((a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(23)/6))
    F = -7776*b**2*d*(a + b*x)**(sympy.S(5)/6)/(935*(c + d*x)**(sympy.S(5)/6)*(-a*d + b*c)**4) - 1296*b*d*(a + b*x)**(sympy.S(5)/6)/(187*(c + d*x)**(sympy.S(11)/6)*(-a*d + b*c)**3) - 108*d*(a + b*x)**(sympy.S(5)/6)/(17*(c + d*x)**(sympy.S(17)/6)*(-a*d + b*c)**2) - 6/((a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(17)/6)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1836():
    f = (c + d*x)**(sympy.S(11)/6)/(a + b*x)**(sympy.S(7)/6)
    F = (c + d*x)**(sympy.S(5)/6)*(6*a*d - 6*b*c)*hyper((sympy.S(-11)/6, sympy.S(-1)/6), (sympy.S(5)/6,), -d*(a + b*x)/(-a*d + b*c))/(b**2*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(5)/6)*(a + b*x)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1837():
    f = (c + d*x)**(sympy.S(5)/6)/(a + b*x)**(sympy.S(7)/6)
    F = -6*(c + d*x)**(sympy.S(5)/6)*hyper((sympy.S(-5)/6, sympy.S(-1)/6), (sympy.S(5)/6,), -d*(a + b*x)/(-a*d + b*c))/(b*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(5)/6)*(a + b*x)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1838():
    f = 1/((a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(1)/6))
    F = -6*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6)*hyper((sympy.S(-1)/6, sympy.S(1)/6), (sympy.S(5)/6,), -d*(a + b*x)/(-a*d + b*c))/(b*(a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1839():
    f = 1/((a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(7)/6))
    F = -6*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6)*hyper((sympy.S(-1)/6, sympy.S(7)/6), (sympy.S(5)/6,), -d*(a + b*x)/(-a*d + b*c))/((a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1840():
    f = 1/((a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(13)/6))
    F = -6*b*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6)*hyper((sympy.S(-1)/6, sympy.S(13)/6), (sympy.S(5)/6,), -d*(a + b*x)/(-a*d + b*c))/((a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1841():
    f = 1/((a + b*x)**(sympy.S(7)/6)*(c + d*x)**(sympy.S(19)/6))
    F = -6*b**2*(b*(c + d*x)/(-a*d + b*c))**(sympy.S(1)/6)*hyper((sympy.S(-1)/6, sympy.S(19)/6), (sympy.S(5)/6,), -d*(a + b*x)/(-a*d + b*c))/((a + b*x)**(sympy.S(1)/6)*(c + d*x)**(sympy.S(1)/6)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1842():
    f = (a + b*x)**m*(a + b*x*(m + 2))
    F = x*(a + b*x)**(m + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1843():
    f = (a + b*x)**m*(c + d*x)**3
    F = d**3*(a + b*x)**(m + 4)/(b**4*(m + 4)) + 3*d**2*(a + b*x)**(m + 3)*(-a*d + b*c)/(b**4*(m + 3)) + 3*d*(a + b*x)**(m + 2)*(-a*d + b*c)**2/(b**4*(m + 2)) + (a + b*x)**(m + 1)*(-a*d + b*c)**3/(b**4*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1844():
    f = (a + b*x)**m*(c + d*x)**2
    F = d**2*(a + b*x)**(m + 3)/(b**3*(m + 3)) + 2*d*(a + b*x)**(m + 2)*(-a*d + b*c)/(b**3*(m + 2)) + (a + b*x)**(m + 1)*(-a*d + b*c)**2/(b**3*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1845():
    f = (a + b*x)**m*(c + d*x)
    F = d*(a + b*x)**(m + 2)/(b**2*(m + 2)) + (a + b*x)**(m + 1)*(-a*d + b*c)/(b**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1846():
    f = (a + b*x)**m/(c + d*x)
    F = (a + b*x)**(m + 1)*hyper((1, m + 1), (m + 2,), -d*(a + b*x)/(-a*d + b*c))/((m + 1)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1847():
    f = (a + b*x)**m/(c + d*x)**2
    F = b*(a + b*x)**(m + 1)*hyper((2, m + 1), (m + 2,), -d*(a + b*x)/(-a*d + b*c))/((m + 1)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1848():
    f = (a + b*x)**m/(c + d*x)**3
    F = b**2*(a + b*x)**(m + 1)*hyper((3, m + 1), (m + 2,), -d*(a + b*x)/(-a*d + b*c))/((m + 1)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1849():
    f = (a + b*x)**3*(c + d*x)**n
    F = b**3*(c + d*x)**(n + 4)/(d**4*(n + 4)) - 3*b**2*(c + d*x)**(n + 3)*(-a*d + b*c)/(d**4*(n + 3)) + 3*b*(c + d*x)**(n + 2)*(-a*d + b*c)**2/(d**4*(n + 2)) - (c + d*x)**(n + 1)*(-a*d + b*c)**3/(d**4*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1850():
    f = (a + b*x)**2*(c + d*x)**n
    F = b**2*(c + d*x)**(n + 3)/(d**3*(n + 3)) - 2*b*(c + d*x)**(n + 2)*(-a*d + b*c)/(d**3*(n + 2)) + (c + d*x)**(n + 1)*(-a*d + b*c)**2/(d**3*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1851():
    f = (a + b*x)*(c + d*x)**n
    F = b*(c + d*x)**(n + 2)/(d**2*(n + 2)) - (c + d*x)**(n + 1)*(-a*d + b*c)/(d**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1852():
    f = (c + d*x)**n
    F = (c + d*x)**(n + 1)/(d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1853():
    f = (c + d*x)**n/(a + b*x)
    F = -(c + d*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b*(c + d*x)/(-a*d + b*c))/((n + 1)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1854():
    f = (c + d*x)**n/(a + b*x)**2
    F = d*(c + d*x)**(n + 1)*hyper((2, n + 1), (n + 2,), b*(c + d*x)/(-a*d + b*c))/((n + 1)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1855():
    f = (c + d*x)**n/(a + b*x)**3
    F = -d**2*(c + d*x)**(n + 1)*hyper((3, n + 1), (n + 2,), b*(c + d*x)/(-a*d + b*c))/((n + 1)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1856():
    f = (a + b*x)**(n - 4)/(c + d*x)**n
    F = -2*d**2*(a + b*x)**(n - 1)*(c + d*x)**(1 - n)/((1 - n)*(2 - n)*(3 - n)*(-a*d + b*c)**3) + 2*d*(a + b*x)**(n - 2)*(c + d*x)**(1 - n)/((2 - n)*(3 - n)*(-a*d + b*c)**2) - (a + b*x)**(n - 3)*(c + d*x)**(1 - n)/((3 - n)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1857():
    f = (a + b*x)**(n - 3)/(c + d*x)**n
    F = d*(a + b*x)**(n - 1)*(c + d*x)**(1 - n)/((1 - n)*(2 - n)*(-a*d + b*c)**2) - (a + b*x)**(n - 2)*(c + d*x)**(1 - n)/((2 - n)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1858():
    f = (a + b*x)**(n - 2)/(c + d*x)**n
    F = -(a + b*x)**(n - 1)*(c + d*x)**(1 - n)/((1 - n)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1859():
    f = (a + b*x)**(n - 1)/(c + d*x)**n
    F = (b*(c + d*x)/(-a*d + b*c))**n*(a + b*x)**n*hyper((n, n), (n + 1,), -d*(a + b*x)/(-a*d + b*c))/(b*n*(c + d*x)**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1860():
    f = (a + b*x)**n/(c + d*x)**n
    F = (b*(c + d*x)/(-a*d + b*c))**n*(a + b*x)**(n + 1)*hyper((n, n + 1), (n + 2,), -d*(a + b*x)/(-a*d + b*c))/(b*(c + d*x)**n*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1861():
    f = (a + b*x)**(n + 1)/(c + d*x)**n
    F = (b*(c + d*x)/(-a*d + b*c))**n*(a + b*x)**(n + 2)*hyper((n, n + 2), (n + 3,), -d*(a + b*x)/(-a*d + b*c))/(b*(c + d*x)**n*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1862():
    f = (a + b*x)**(n + 2)/(c + d*x)**n
    F = (b*(c + d*x)/(-a*d + b*c))**n*(a + b*x)**(n + 3)*hyper((n, n + 3), (n + 4,), -d*(a + b*x)/(-a*d + b*c))/(b*(c + d*x)**n*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1863():
    f = (c + d*x)**n/(a + b*x)**n
    F = (-d*(a + b*x)/(-a*d + b*c))**n*(c + d*x)**(n + 1)*hyper((n, n + 1), (n + 2,), b*(c + d*x)/(-a*d + b*c))/(d*(a + b*x)**n*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1864():
    f = (a + b*x)**(-n - 1)*(c + d*x)**n
    F = -(c + d*x)**n*hyper((-n, -n), (1 - n,), -d*(a + b*x)/(-a*d + b*c))/(b*n*(b*(c + d*x)/(-a*d + b*c))**n*(a + b*x)**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1865():
    f = (a + b*x)**(-n - 2)*(c + d*x)**n
    F = -(a + b*x)**(-n - 1)*(c + d*x)**(n + 1)/((n + 1)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1866():
    f = (a + b*x)**(-n - 3)*(c + d*x)**n
    F = d*(a + b*x)**(-n - 1)*(c + d*x)**(n + 1)/((n + 1)*(n + 2)*(-a*d + b*c)**2) - (a + b*x)**(-n - 2)*(c + d*x)**(n + 1)/((n + 2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1867():
    f = (a + b*x)**(-n - 4)*(c + d*x)**n
    F = -2*d**2*(a + b*x)**(-n - 1)*(c + d*x)**(n + 1)/((n + 1)*(n + 2)*(n + 3)*(-a*d + b*c)**3) + 2*d*(a + b*x)**(-n - 2)*(c + d*x)**(n + 1)/((n + 2)*(n + 3)*(-a*d + b*c)**2) - (a + b*x)**(-n - 3)*(c + d*x)**(n + 1)/((n + 3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1868():
    f = (a + b*x)**(-n - 5)*(c + d*x)**n
    F = 6*d**3*(a + b*x)**(-n - 1)*(c + d*x)**(n + 1)/((n + 1)*(n + 2)*(n + 3)*(n + 4)*(-a*d + b*c)**4) - 6*d**2*(a + b*x)**(-n - 2)*(c + d*x)**(n + 1)/((n + 2)*(n + 3)*(n + 4)*(-a*d + b*c)**3) + 3*d*(a + b*x)**(-n - 3)*(c + d*x)**(n + 1)/((n + 3)*(n + 4)*(-a*d + b*c)**2) - (a + b*x)**(-n - 4)*(c + d*x)**(n + 1)/((n + 4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1869():
    f = (a + b*x)**n/(c + d*x)**n
    F = (b*(c + d*x)/(-a*d + b*c))**n*(a + b*x)**(n + 1)*hyper((n, n + 1), (n + 2,), -d*(a + b*x)/(-a*d + b*c))/(b*(c + d*x)**n*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1870():
    f = (a + b*x)**n*(c + d*x)**(-n - 1)
    F = -(a + b*x)**n*hyper((-n, -n), (1 - n,), b*(c + d*x)/(-a*d + b*c))/(d*n*(-d*(a + b*x)/(-a*d + b*c))**n*(c + d*x)**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1871():
    f = (a + b*x)**n*(c + d*x)**(-n - 2)
    F = (a + b*x)**(n + 1)*(c + d*x)**(-n - 1)/((n + 1)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1872():
    f = (a + b*x)**n*(c + d*x)**(-n - 3)
    F = b*(a + b*x)**(n + 1)*(c + d*x)**(-n - 1)/((n + 1)*(n + 2)*(-a*d + b*c)**2) + (a + b*x)**(n + 1)*(c + d*x)**(-n - 2)/((n + 2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1873():
    f = (a + b*x)**n*(c + d*x)**(-n - 4)
    F = 2*b**2*(a + b*x)**(n + 1)*(c + d*x)**(-n - 1)/((n + 1)*(n + 2)*(n + 3)*(-a*d + b*c)**3) + 2*b*(a + b*x)**(n + 1)*(c + d*x)**(-n - 2)/((n + 2)*(n + 3)*(-a*d + b*c)**2) + (a + b*x)**(n + 1)*(c + d*x)**(-n - 3)/((n + 3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1874():
    f = (a + b*x)**n*(c + d*x)**(-n - 5)
    F = 6*b**3*(a + b*x)**(n + 1)*(c + d*x)**(-n - 1)/((n + 1)*(n + 2)*(n + 3)*(n + 4)*(-a*d + b*c)**4) + 6*b**2*(a + b*x)**(n + 1)*(c + d*x)**(-n - 2)/((n + 2)*(n + 3)*(n + 4)*(-a*d + b*c)**3) + 3*b*(a + b*x)**(n + 1)*(c + d*x)**(-n - 3)/((n + 3)*(n + 4)*(-a*d + b*c)**2) + (a + b*x)**(n + 1)*(c + d*x)**(-n - 4)/((n + 4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1875():
    f = (a + b*x)**(n - 2)*(c + d*x)**(1 - n)
    F = -(b*(c + d*x)/(-a*d + b*c))**n*(a + b*x)**(n - 1)*(-a*d + b*c)*hyper((n - 1, n - 1), (n,), -d*(a + b*x)/(-a*d + b*c))/(b**2*(1 - n)*(c + d*x)**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1876():
    f = (a + b*x)**(n + 1)*(c + d*x)**(-n - 1)
    F = (a + b*x)**n*(-a*d + b*c)*hyper((-n, -n - 1), (1 - n,), b*(c + d*x)/(-a*d + b*c))/(d**2*n*(-d*(a + b*x)/(-a*d + b*c))**n*(c + d*x)**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1877():
    f = (a + b*x)**m/(c + d*x)
    F = (a + b*x)**(m + 1)*hyper((1, m + 1), (m + 2,), -d*(a + b*x)/(-a*d + b*c))/((m + 1)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1878():
    f = 1/((a + b*x)**2*(c + d*x))
    F = -d*log(a + b*x)/(-a*d + b*c)**2 + d*log(c + d*x)/(-a*d + b*c)**2 - 1/((a + b*x)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1879():
    f = (a + b*x)**m*(a*c*(m + 1) + b*c*x*(m + 2))**(-m - 3)
    F = -(a + b*x)**(m + 1)*(a*c*(m + 1) + b*c*x*(m + 2))**(-m - 2)/(a*b*c*(m + 2)) + (a + b*x)**(m + 1)*(a*c*(m + 1) + b*c*x*(m + 2))**(-m - 1)/(a**2*b*c**2*(m + 1)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1880():
    f = (a + b*x)**(-b*c/(-a*d + b*c) - 1)*(c + d*x)**(a*d/(-a*d + b*c) - 1)
    F = -(c + d*x)**(a*d/(-a*d + b*c))/(b*c*(a + b*x)**(b*c/(-a*d + b*c))) + (c + d*x)**(a*d/(-a*d + b*c))/(a*b*c*(a + b*x)**(a*d/(-a*d + b*c)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1881():
    f = (a + b*x)**((a*d - 2*b*c)/(-a*d + b*c))*(c + d*x)**((-2*a*d + b*c)/(a*d - b*c))
    F = -(c + d*x)**(a*d/(-a*d + b*c))/(b*c*(a + b*x)**(b*c/(-a*d + b*c))) + (c + d*x)**(a*d/(-a*d + b*c))/(a*b*c*(a + b*x)**(a*d/(-a*d + b*c)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1882():
    f = (1 - x)**n/sqrt(x + 1)
    F = 2**(n + 1)*sqrt(x + 1)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), x/2 + sympy.S.Half)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1883():
    f = (x + 1)**n/sqrt(1 - x)
    F = -2**(n + 1)*sqrt(1 - x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), sympy.S.Half - x/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1884():
    f = (1 - x)**n*(x + 1)**(sympy.S(7)/3)
    F = 3*2**(n - 1)*(x + 1)**(sympy.S(10)/3)*hyper((sympy.S(10)/3, -n), (sympy.S(13)/3,), x/2 + sympy.S.Half)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1885():
    f = (1 - x)**(sympy.S(7)/3)*(x + 1)**n
    F = -3*2**(n - 1)*(1 - x)**(sympy.S(10)/3)*hyper((sympy.S(10)/3, -n), (sympy.S(13)/3,), sympy.S.Half - x/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1886():
    f = (3*x + 2)**m/(2*x + 1)**m
    F = 2**(-m - 1)*(2*x + 1)**(1 - m)*hyper((-m, 1 - m), (2 - m,), -6*x - 3)/(1 - m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1887():
    f = (d*(a + b*x)/(a*d - b*c))**m*(c + d*x)**n
    F = (c + d*x)**(n + 1)*hyper((-m, n + 1), (n + 2,), b*(c + d*x)/(-a*d + b*c))/(d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1888():
    f = a + b*x + c*x**2 + d*x**3
    F = a*x + b*x**2/2 + c*x**3/3 + d*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1889():
    f = x**4 - x**3
    F = x**5/5 - x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1890():
    f = x**5 - 1
    F = x**6/6 - x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1891():
    f = 4*x + 7
    F = 2*x**2 + 7*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1892():
    f = pi*x**3 + 4*x
    F = pi*x**4/4 + 2*x**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1893():
    f = 5*x**2 + 2*x
    F = 5*x**3/3 + x**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1894():
    f = x**3/3 + x**2/2
    F = x**4/12 + x**3/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1895():
    f = 2*x**2 - 5*x + 3
    F = 2*x**3/3 - 5*x**2/2 + 3*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1896():
    f = x**3 + x**2 - 2*x
    F = x**4/4 + x**3/3 - x**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1897():
    f = -3*x**5 - x**2 + 1
    F = -x**6/2 - x**3/3 + x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1898():
    f = 4*x**3 + 3*x**2 + 2*x + 5
    F = x**4 + x**3 + x**2 + 5*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1899():
    f = a + b/x + c/x**2 + d/x**3
    F = a*x + b*log(x) - c/x - d/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1900():
    f = x**5 + x + x**(-5)
    F = x**6/6 + x**2/2 - 1/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1901():
    f = 1/x + x**(-2) + x**(-3)
    F = log(x) - 1/x - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1902():
    f = 3/x - 2/x**2
    F = 3*log(x) + 2/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1903():
    f = x**6 - 1/(7*x**6)
    F = x**7/7 + 1/(35*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1904():
    f = x + 1 + 1/x
    F = x**2/2 + x + log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1905():
    f = 4/x**2 - 3/x**3
    F = -4/x + 3/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1906():
    f = x**2 + 2*x + 1/x
    F = x**3/3 + x**2 + log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1907():
    f = x**(sympy.S(5)/6) - x**3
    F = 6*x**(sympy.S(11)/6)/11 - x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1908():
    f = x**(sympy.S(1)/33) + 33
    F = 33*x**(sympy.S(34)/33)/34 + 33*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1909():
    f = 2*sqrt(x) + 1/(2*sqrt(x))
    F = 4*x**(sympy.S(3)/2)/3 + sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1910():
    f = 6*sqrt(x) + 10/x - 1/x**2
    F = 4*x**(sympy.S(3)/2) + 10*log(x) + 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1911():
    f = x**(sympy.S(3)/2) + x**(sympy.S(-3)/2)
    F = 2*x**(sympy.S(5)/2)/5 - 2/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1912():
    f = 7*x**(sympy.S(5)/2) - 5*x**(sympy.S(3)/2)
    F = 2*x**(sympy.S(7)/2) - 2*x**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1913():
    f = sqrt(x) - x/2 + 2/sqrt(x)
    F = 2*x**(sympy.S(3)/2)/3 + 4*sqrt(x) - x**2/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_2_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1914():
    f = x**(sympy.S(3)/2) + sqrt(x)/5 - 2/x
    F = 2*x**(sympy.S(5)/2)/5 + 2*x**(sympy.S(3)/2)/15 - 2*log(x)
    assert integrate(f, x) == F

