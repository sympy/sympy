"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.1 Binomial products/1.1.2 Quadratic/1.1.2.2 (c x)^m (a+b x^2)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, m, p = symbols('a b c d m p')

def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1():
    f = x**4*(a + b*x**2)
    F = a*x**5/5 + b*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_2():
    f = x**3*(a + b*x**2)
    F = a*x**4/4 + b*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_3():
    f = x**2*(a + b*x**2)
    F = a*x**3/3 + b*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_4():
    f = x*(a + b*x**2)
    F = a*x**2/2 + b*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_5():
    f = a + b*x**2
    F = a*x + b*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_6():
    f = (a + b*x**2)/x
    F = a*log(x) + b*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_7():
    f = (a + b*x**2)/x**2
    F = -a/x + b*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_8():
    f = (a + b*x**2)/x**3
    F = -a/(2*x**2) + b*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_9():
    f = (a + b*x**2)/x**4
    F = -a/(3*x**3) - b/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_10():
    f = (a + b*x**2)/x**5
    F = -a/(4*x**4) - b/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_11():
    f = (a + b*x**2)/x**6
    F = -a/(5*x**5) - b/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_12():
    f = (a + b*x**2)/x**7
    F = -a/(6*x**6) - b/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_13():
    f = x**5*(a + b*x**2)**2
    F = a**2*x**6/6 + a*b*x**8/4 + b**2*x**10/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_14():
    f = x**4*(a + b*x**2)**2
    F = a**2*x**5/5 + 2*a*b*x**7/7 + b**2*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_15():
    f = x**3*(a + b*x**2)**2
    F = a**2*x**4/4 + a*b*x**6/3 + b**2*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_16():
    f = x**2*(a + b*x**2)**2
    F = a**2*x**3/3 + 2*a*b*x**5/5 + b**2*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_17():
    f = x*(a + b*x**2)**2
    F = (a + b*x**2)**3/(6*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_18():
    f = (a + b*x**2)**2
    F = a**2*x + 2*a*b*x**3/3 + b**2*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_19():
    f = (a + b*x**2)**2/x
    F = a**2*log(x) + a*b*x**2 + b**2*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_20():
    f = (a + b*x**2)**2/x**2
    F = -a**2/x + 2*a*b*x + b**2*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_21():
    f = (a + b*x**2)**2/x**3
    F = -a**2/(2*x**2) + 2*a*b*log(x) + b**2*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_22():
    f = (a + b*x**2)**2/x**4
    F = -a**2/(3*x**3) - 2*a*b/x + b**2*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_23():
    f = (a + b*x**2)**2/x**5
    F = -a**2/(4*x**4) - a*b/x**2 + b**2*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_24():
    f = (a + b*x**2)**2/x**6
    F = -a**2/(5*x**5) - 2*a*b/(3*x**3) - b**2/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_25():
    f = (a + b*x**2)**2/x**7
    F = -(a + b*x**2)**3/(6*a*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_26():
    f = (a + b*x**2)**2/x**8
    F = -a**2/(7*x**7) - 2*a*b/(5*x**5) - b**2/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_27():
    f = (a + b*x**2)**2/x**9
    F = -a**2/(8*x**8) - a*b/(3*x**6) - b**2/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_28():
    f = (a + b*x**2)**2/x**10
    F = -a**2/(9*x**9) - 2*a*b/(7*x**7) - b**2/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_29():
    f = x**9*(a + b*x**2)**3
    F = a**3*x**10/10 + a**2*b*x**12/4 + 3*a*b**2*x**14/14 + b**3*x**16/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_30():
    f = x**7*(a + b*x**2)**3
    F = a**3*x**8/8 + 3*a**2*b*x**10/10 + a*b**2*x**12/4 + b**3*x**14/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_31():
    f = x**5*(a + b*x**2)**3
    F = a**3*x**6/6 + 3*a**2*b*x**8/8 + 3*a*b**2*x**10/10 + b**3*x**12/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_32():
    f = x**3*(a + b*x**2)**3
    F = -a*(a + b*x**2)**4/(8*b**2) + (a + b*x**2)**5/(10*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_33():
    f = x*(a + b*x**2)**3
    F = (a + b*x**2)**4/(8*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_34():
    f = (a + b*x**2)**3/x
    F = a**3*log(x) + 3*a**2*b*x**2/2 + 3*a*b**2*x**4/4 + b**3*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_35():
    f = (a + b*x**2)**3/x**3
    F = -a**3/(2*x**2) + 3*a**2*b*log(x) + 3*a*b**2*x**2/2 + b**3*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_36():
    f = (a + b*x**2)**3/x**5
    F = -a**3/(4*x**4) - 3*a**2*b/(2*x**2) + 3*a*b**2*log(x) + b**3*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_37():
    f = (a + b*x**2)**3/x**7
    F = -a**3/(6*x**6) - 3*a**2*b/(4*x**4) - 3*a*b**2/(2*x**2) + b**3*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_38():
    f = (a + b*x**2)**3/x**9
    F = -(a + b*x**2)**4/(8*a*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_39():
    f = (a + b*x**2)**3/x**11
    F = -(a + b*x**2)**4/(10*a*x**10) + b*(a + b*x**2)**4/(40*a**2*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_40():
    f = (a + b*x**2)**3/x**13
    F = -a**3/(12*x**12) - 3*a**2*b/(10*x**10) - 3*a*b**2/(8*x**8) - b**3/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_41():
    f = (a + b*x**2)**3/x**15
    F = -a**3/(14*x**14) - a**2*b/(4*x**12) - 3*a*b**2/(10*x**10) - b**3/(8*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_42():
    f = x**6*(a + b*x**2)**3
    F = a**3*x**7/7 + a**2*b*x**9/3 + 3*a*b**2*x**11/11 + b**3*x**13/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_43():
    f = x**4*(a + b*x**2)**3
    F = a**3*x**5/5 + 3*a**2*b*x**7/7 + a*b**2*x**9/3 + b**3*x**11/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_44():
    f = x**2*(a + b*x**2)**3
    F = a**3*x**3/3 + 3*a**2*b*x**5/5 + 3*a*b**2*x**7/7 + b**3*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_45():
    f = (a + b*x**2)**3
    F = a**3*x + a**2*b*x**3 + 3*a*b**2*x**5/5 + b**3*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_46():
    f = (a + b*x**2)**3/x**2
    F = -a**3/x + 3*a**2*b*x + a*b**2*x**3 + b**3*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_47():
    f = (a + b*x**2)**3/x**4
    F = -a**3/(3*x**3) - 3*a**2*b/x + 3*a*b**2*x + b**3*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_48():
    f = (a + b*x**2)**3/x**6
    F = -a**3/(5*x**5) - a**2*b/x**3 - 3*a*b**2/x + b**3*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_49():
    f = (a + b*x**2)**3/x**8
    F = -a**3/(7*x**7) - 3*a**2*b/(5*x**5) - a*b**2/x**3 - b**3/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_50():
    f = (a + b*x**2)**3/x**10
    F = -a**3/(9*x**9) - 3*a**2*b/(7*x**7) - 3*a*b**2/(5*x**5) - b**3/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_51():
    f = (a + b*x**2)**3/x**12
    F = -a**3/(11*x**11) - a**2*b/(3*x**9) - 3*a*b**2/(7*x**7) - b**3/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_52():
    f = x**13*(a + b*x**2)**5
    F = a**5*x**14/14 + 5*a**4*b*x**16/16 + 5*a**3*b**2*x**18/9 + a**2*b**3*x**20/2 + 5*a*b**4*x**22/22 + b**5*x**24/24
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_53():
    f = x**11*(a + b*x**2)**5
    F = a**5*x**12/12 + 5*a**4*b*x**14/14 + 5*a**3*b**2*x**16/8 + 5*a**2*b**3*x**18/9 + a*b**4*x**20/4 + b**5*x**22/22
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_54():
    f = x**9*(a + b*x**2)**5
    F = a**5*x**10/10 + 5*a**4*b*x**12/12 + 5*a**3*b**2*x**14/7 + 5*a**2*b**3*x**16/8 + 5*a*b**4*x**18/18 + b**5*x**20/20
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_55():
    f = x**7*(a + b*x**2)**5
    F = -a**3*(a + b*x**2)**6/(12*b**4) + 3*a**2*(a + b*x**2)**7/(14*b**4) - 3*a*(a + b*x**2)**8/(16*b**4) + (a + b*x**2)**9/(18*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_56():
    f = x**5*(a + b*x**2)**5
    F = a**2*(a + b*x**2)**6/(12*b**3) - a*(a + b*x**2)**7/(7*b**3) + (a + b*x**2)**8/(16*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_57():
    f = x**3*(a + b*x**2)**5
    F = -a*(a + b*x**2)**6/(12*b**2) + (a + b*x**2)**7/(14*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_58():
    f = x*(a + b*x**2)**5
    F = (a + b*x**2)**6/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_59():
    f = (a + b*x**2)**5/x
    F = a**5*log(x) + 5*a**4*b*x**2/2 + 5*a**3*b**2*x**4/2 + 5*a**2*b**3*x**6/3 + 5*a*b**4*x**8/8 + b**5*x**10/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_60():
    f = (a + b*x**2)**5/x**3
    F = -a**5/(2*x**2) + 5*a**4*b*log(x) + 5*a**3*b**2*x**2 + 5*a**2*b**3*x**4/2 + 5*a*b**4*x**6/6 + b**5*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_61():
    f = (a + b*x**2)**5/x**5
    F = -a**5/(4*x**4) - 5*a**4*b/(2*x**2) + 10*a**3*b**2*log(x) + 5*a**2*b**3*x**2 + 5*a*b**4*x**4/4 + b**5*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_62():
    f = (a + b*x**2)**5/x**7
    F = -a**5/(6*x**6) - 5*a**4*b/(4*x**4) - 5*a**3*b**2/x**2 + 10*a**2*b**3*log(x) + 5*a*b**4*x**2/2 + b**5*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_63():
    f = (a + b*x**2)**5/x**9
    F = -a**5/(8*x**8) - 5*a**4*b/(6*x**6) - 5*a**3*b**2/(2*x**4) - 5*a**2*b**3/x**2 + 5*a*b**4*log(x) + b**5*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_64():
    f = (a + b*x**2)**5/x**11
    F = -a**5/(10*x**10) - 5*a**4*b/(8*x**8) - 5*a**3*b**2/(3*x**6) - 5*a**2*b**3/(2*x**4) - 5*a*b**4/(2*x**2) + b**5*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_65():
    f = (a + b*x**2)**5/x**13
    F = -(a + b*x**2)**6/(12*a*x**12)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_66():
    f = (a + b*x**2)**5/x**15
    F = -(a + b*x**2)**6/(14*a*x**14) + b*(a + b*x**2)**6/(84*a**2*x**12)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_67():
    f = (a + b*x**2)**5/x**17
    F = -(a + b*x**2)**6/(16*a*x**16) + b*(a + b*x**2)**6/(56*a**2*x**14) - b**2*(a + b*x**2)**6/(336*a**3*x**12)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_68():
    f = (a + b*x**2)**5/x**19
    F = -a**5/(18*x**18) - 5*a**4*b/(16*x**16) - 5*a**3*b**2/(7*x**14) - 5*a**2*b**3/(6*x**12) - a*b**4/(2*x**10) - b**5/(8*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_69():
    f = (a + b*x**2)**5/x**21
    F = -a**5/(20*x**20) - 5*a**4*b/(18*x**18) - 5*a**3*b**2/(8*x**16) - 5*a**2*b**3/(7*x**14) - 5*a*b**4/(12*x**12) - b**5/(10*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_70():
    f = x**8*(a + b*x**2)**5
    F = a**5*x**9/9 + 5*a**4*b*x**11/11 + 10*a**3*b**2*x**13/13 + 2*a**2*b**3*x**15/3 + 5*a*b**4*x**17/17 + b**5*x**19/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_71():
    f = x**6*(a + b*x**2)**5
    F = a**5*x**7/7 + 5*a**4*b*x**9/9 + 10*a**3*b**2*x**11/11 + 10*a**2*b**3*x**13/13 + a*b**4*x**15/3 + b**5*x**17/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_72():
    f = x**4*(a + b*x**2)**5
    F = a**5*x**5/5 + 5*a**4*b*x**7/7 + 10*a**3*b**2*x**9/9 + 10*a**2*b**3*x**11/11 + 5*a*b**4*x**13/13 + b**5*x**15/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_73():
    f = x**2*(a + b*x**2)**5
    F = a**5*x**3/3 + a**4*b*x**5 + 10*a**3*b**2*x**7/7 + 10*a**2*b**3*x**9/9 + 5*a*b**4*x**11/11 + b**5*x**13/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_74():
    f = (a + b*x**2)**5
    F = a**5*x + 5*a**4*b*x**3/3 + 2*a**3*b**2*x**5 + 10*a**2*b**3*x**7/7 + 5*a*b**4*x**9/9 + b**5*x**11/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_75():
    f = (a + b*x**2)**5/x**2
    F = -a**5/x + 5*a**4*b*x + 10*a**3*b**2*x**3/3 + 2*a**2*b**3*x**5 + 5*a*b**4*x**7/7 + b**5*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_76():
    f = (a + b*x**2)**5/x**4
    F = -a**5/(3*x**3) - 5*a**4*b/x + 10*a**3*b**2*x + 10*a**2*b**3*x**3/3 + a*b**4*x**5 + b**5*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_77():
    f = (a + b*x**2)**5/x**6
    F = -a**5/(5*x**5) - 5*a**4*b/(3*x**3) - 10*a**3*b**2/x + 10*a**2*b**3*x + 5*a*b**4*x**3/3 + b**5*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_78():
    f = (a + b*x**2)**5/x**8
    F = -a**5/(7*x**7) - a**4*b/x**5 - 10*a**3*b**2/(3*x**3) - 10*a**2*b**3/x + 5*a*b**4*x + b**5*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_79():
    f = (a + b*x**2)**5/x**10
    F = -a**5/(9*x**9) - 5*a**4*b/(7*x**7) - 2*a**3*b**2/x**5 - 10*a**2*b**3/(3*x**3) - 5*a*b**4/x + b**5*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_80():
    f = (a + b*x**2)**5/x**12
    F = -a**5/(11*x**11) - 5*a**4*b/(9*x**9) - 10*a**3*b**2/(7*x**7) - 2*a**2*b**3/x**5 - 5*a*b**4/(3*x**3) - b**5/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_81():
    f = (a + b*x**2)**5/x**14
    F = -a**5/(13*x**13) - 5*a**4*b/(11*x**11) - 10*a**3*b**2/(9*x**9) - 10*a**2*b**3/(7*x**7) - a*b**4/x**5 - b**5/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_82():
    f = (a + b*x**2)**5/x**16
    F = -a**5/(15*x**15) - 5*a**4*b/(13*x**13) - 10*a**3*b**2/(11*x**11) - 10*a**2*b**3/(9*x**9) - 5*a*b**4/(7*x**7) - b**5/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_83():
    f = (a + b*x**2)**5/x**18
    F = -a**5/(17*x**17) - a**4*b/(3*x**15) - 10*a**3*b**2/(13*x**13) - 10*a**2*b**3/(11*x**11) - 5*a*b**4/(9*x**9) - b**5/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_84():
    f = (a + b*x**2)**5/x**20
    F = -a**5/(19*x**19) - 5*a**4*b/(17*x**17) - 2*a**3*b**2/(3*x**15) - 10*a**2*b**3/(13*x**13) - 5*a*b**4/(11*x**11) - b**5/(9*x**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_85():
    f = x**13*(a + b*x**2)**8
    F = a**6*(a + b*x**2)**9/(18*b**7) - 3*a**5*(a + b*x**2)**10/(10*b**7) + 15*a**4*(a + b*x**2)**11/(22*b**7) - 5*a**3*(a + b*x**2)**12/(6*b**7) + 15*a**2*(a + b*x**2)**13/(26*b**7) - 3*a*(a + b*x**2)**14/(14*b**7) + (a + b*x**2)**15/(30*b**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_86():
    f = x**11*(a + b*x**2)**8
    F = -a**5*(a + b*x**2)**9/(18*b**6) + a**4*(a + b*x**2)**10/(4*b**6) - 5*a**3*(a + b*x**2)**11/(11*b**6) + 5*a**2*(a + b*x**2)**12/(12*b**6) - 5*a*(a + b*x**2)**13/(26*b**6) + (a + b*x**2)**14/(28*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_87():
    f = x**9*(a + b*x**2)**8
    F = a**4*(a + b*x**2)**9/(18*b**5) - a**3*(a + b*x**2)**10/(5*b**5) + 3*a**2*(a + b*x**2)**11/(11*b**5) - a*(a + b*x**2)**12/(6*b**5) + (a + b*x**2)**13/(26*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_88():
    f = x**7*(a + b*x**2)**8
    F = -a**3*(a + b*x**2)**9/(18*b**4) + 3*a**2*(a + b*x**2)**10/(20*b**4) - 3*a*(a + b*x**2)**11/(22*b**4) + (a + b*x**2)**12/(24*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_89():
    f = x**5*(a + b*x**2)**8
    F = a**2*(a + b*x**2)**9/(18*b**3) - a*(a + b*x**2)**10/(10*b**3) + (a + b*x**2)**11/(22*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_90():
    f = x**3*(a + b*x**2)**8
    F = -a*(a + b*x**2)**9/(18*b**2) + (a + b*x**2)**10/(20*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_91():
    f = x*(a + b*x**2)**8
    F = (a + b*x**2)**9/(18*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_92():
    f = (a + b*x**2)**8/x
    F = a**8*log(x) + 4*a**7*b*x**2 + 7*a**6*b**2*x**4 + 28*a**5*b**3*x**6/3 + 35*a**4*b**4*x**8/4 + 28*a**3*b**5*x**10/5 + 7*a**2*b**6*x**12/3 + 4*a*b**7*x**14/7 + b**8*x**16/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_93():
    f = (a + b*x**2)**8/x**3
    F = -a**8/(2*x**2) + 8*a**7*b*log(x) + 14*a**6*b**2*x**2 + 14*a**5*b**3*x**4 + 35*a**4*b**4*x**6/3 + 7*a**3*b**5*x**8 + 14*a**2*b**6*x**10/5 + 2*a*b**7*x**12/3 + b**8*x**14/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_94():
    f = (a + b*x**2)**8/x**5
    F = -a**8/(4*x**4) - 4*a**7*b/x**2 + 28*a**6*b**2*log(x) + 28*a**5*b**3*x**2 + 35*a**4*b**4*x**4/2 + 28*a**3*b**5*x**6/3 + 7*a**2*b**6*x**8/2 + 4*a*b**7*x**10/5 + b**8*x**12/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_95():
    f = (a + b*x**2)**8/x**7
    F = -a**8/(6*x**6) - 2*a**7*b/x**4 - 14*a**6*b**2/x**2 + 56*a**5*b**3*log(x) + 35*a**4*b**4*x**2 + 14*a**3*b**5*x**4 + 14*a**2*b**6*x**6/3 + a*b**7*x**8 + b**8*x**10/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_96():
    f = (a + b*x**2)**8/x**9
    F = -a**8/(8*x**8) - 4*a**7*b/(3*x**6) - 7*a**6*b**2/x**4 - 28*a**5*b**3/x**2 + 70*a**4*b**4*log(x) + 28*a**3*b**5*x**2 + 7*a**2*b**6*x**4 + 4*a*b**7*x**6/3 + b**8*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_97():
    f = (a + b*x**2)**8/x**11
    F = -a**8/(10*x**10) - a**7*b/x**8 - 14*a**6*b**2/(3*x**6) - 14*a**5*b**3/x**4 - 35*a**4*b**4/x**2 + 56*a**3*b**5*log(x) + 14*a**2*b**6*x**2 + 2*a*b**7*x**4 + b**8*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_98():
    f = (a + b*x**2)**8/x**13
    F = -a**8/(12*x**12) - 4*a**7*b/(5*x**10) - 7*a**6*b**2/(2*x**8) - 28*a**5*b**3/(3*x**6) - 35*a**4*b**4/(2*x**4) - 28*a**3*b**5/x**2 + 28*a**2*b**6*log(x) + 4*a*b**7*x**2 + b**8*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_99():
    f = (a + b*x**2)**8/x**15
    F = -a**8/(14*x**14) - 2*a**7*b/(3*x**12) - 14*a**6*b**2/(5*x**10) - 7*a**5*b**3/x**8 - 35*a**4*b**4/(3*x**6) - 14*a**3*b**5/x**4 - 14*a**2*b**6/x**2 + 8*a*b**7*log(x) + b**8*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_100():
    f = (a + b*x**2)**8/x**17
    F = -a**8/(16*x**16) - 4*a**7*b/(7*x**14) - 7*a**6*b**2/(3*x**12) - 28*a**5*b**3/(5*x**10) - 35*a**4*b**4/(4*x**8) - 28*a**3*b**5/(3*x**6) - 7*a**2*b**6/x**4 - 4*a*b**7/x**2 + b**8*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_101():
    f = (a + b*x**2)**8/x**19
    F = -(a + b*x**2)**9/(18*a*x**18)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_102():
    f = (a + b*x**2)**8/x**21
    F = -(a + b*x**2)**9/(20*a*x**20) + b*(a + b*x**2)**9/(180*a**2*x**18)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_103():
    f = (a + b*x**2)**8/x**23
    F = -(a + b*x**2)**9/(22*a*x**22) + b*(a + b*x**2)**9/(110*a**2*x**20) - b**2*(a + b*x**2)**9/(990*a**3*x**18)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_104():
    f = (a + b*x**2)**8/x**25
    F = -(a + b*x**2)**9/(24*a*x**24) + b*(a + b*x**2)**9/(88*a**2*x**22) - b**2*(a + b*x**2)**9/(440*a**3*x**20) + b**3*(a + b*x**2)**9/(3960*a**4*x**18)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_105():
    f = (a + b*x**2)**8/x**27
    F = -(a + b*x**2)**9/(26*a*x**26) + b*(a + b*x**2)**9/(78*a**2*x**24) - b**2*(a + b*x**2)**9/(286*a**3*x**22) + b**3*(a + b*x**2)**9/(1430*a**4*x**20) - b**4*(a + b*x**2)**9/(12870*a**5*x**18)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_106():
    f = (a + b*x**2)**8/x**29
    F = -a**8/(28*x**28) - 4*a**7*b/(13*x**26) - 7*a**6*b**2/(6*x**24) - 28*a**5*b**3/(11*x**22) - 7*a**4*b**4/(2*x**20) - 28*a**3*b**5/(9*x**18) - 7*a**2*b**6/(4*x**16) - 4*a*b**7/(7*x**14) - b**8/(12*x**12)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_107():
    f = (a + b*x**2)**8/x**31
    F = -a**8/(30*x**30) - 2*a**7*b/(7*x**28) - 14*a**6*b**2/(13*x**26) - 7*a**5*b**3/(3*x**24) - 35*a**4*b**4/(11*x**22) - 14*a**3*b**5/(5*x**20) - 14*a**2*b**6/(9*x**18) - a*b**7/(2*x**16) - b**8/(14*x**14)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_108():
    f = (a + b*x**2)**8/x**33
    F = -a**8/(32*x**32) - 4*a**7*b/(15*x**30) - a**6*b**2/x**28 - 28*a**5*b**3/(13*x**26) - 35*a**4*b**4/(12*x**24) - 28*a**3*b**5/(11*x**22) - 7*a**2*b**6/(5*x**20) - 4*a*b**7/(9*x**18) - b**8/(16*x**16)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_109():
    f = x**8*(a + b*x**2)**8
    F = a**8*x**9/9 + 8*a**7*b*x**11/11 + 28*a**6*b**2*x**13/13 + 56*a**5*b**3*x**15/15 + 70*a**4*b**4*x**17/17 + 56*a**3*b**5*x**19/19 + 4*a**2*b**6*x**21/3 + 8*a*b**7*x**23/23 + b**8*x**25/25
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_110():
    f = x**6*(a + b*x**2)**8
    F = a**8*x**7/7 + 8*a**7*b*x**9/9 + 28*a**6*b**2*x**11/11 + 56*a**5*b**3*x**13/13 + 14*a**4*b**4*x**15/3 + 56*a**3*b**5*x**17/17 + 28*a**2*b**6*x**19/19 + 8*a*b**7*x**21/21 + b**8*x**23/23
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_111():
    f = x**4*(a + b*x**2)**8
    F = a**8*x**5/5 + 8*a**7*b*x**7/7 + 28*a**6*b**2*x**9/9 + 56*a**5*b**3*x**11/11 + 70*a**4*b**4*x**13/13 + 56*a**3*b**5*x**15/15 + 28*a**2*b**6*x**17/17 + 8*a*b**7*x**19/19 + b**8*x**21/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_112():
    f = x**2*(a + b*x**2)**8
    F = a**8*x**3/3 + 8*a**7*b*x**5/5 + 4*a**6*b**2*x**7 + 56*a**5*b**3*x**9/9 + 70*a**4*b**4*x**11/11 + 56*a**3*b**5*x**13/13 + 28*a**2*b**6*x**15/15 + 8*a*b**7*x**17/17 + b**8*x**19/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_113():
    f = (a + b*x**2)**8
    F = a**8*x + 8*a**7*b*x**3/3 + 28*a**6*b**2*x**5/5 + 8*a**5*b**3*x**7 + 70*a**4*b**4*x**9/9 + 56*a**3*b**5*x**11/11 + 28*a**2*b**6*x**13/13 + 8*a*b**7*x**15/15 + b**8*x**17/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_114():
    f = (a + b*x**2)**8/x**2
    F = -a**8/x + 8*a**7*b*x + 28*a**6*b**2*x**3/3 + 56*a**5*b**3*x**5/5 + 10*a**4*b**4*x**7 + 56*a**3*b**5*x**9/9 + 28*a**2*b**6*x**11/11 + 8*a*b**7*x**13/13 + b**8*x**15/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_115():
    f = (a + b*x**2)**8/x**4
    F = -a**8/(3*x**3) - 8*a**7*b/x + 28*a**6*b**2*x + 56*a**5*b**3*x**3/3 + 14*a**4*b**4*x**5 + 8*a**3*b**5*x**7 + 28*a**2*b**6*x**9/9 + 8*a*b**7*x**11/11 + b**8*x**13/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_116():
    f = (a + b*x**2)**8/x**6
    F = -a**8/(5*x**5) - 8*a**7*b/(3*x**3) - 28*a**6*b**2/x + 56*a**5*b**3*x + 70*a**4*b**4*x**3/3 + 56*a**3*b**5*x**5/5 + 4*a**2*b**6*x**7 + 8*a*b**7*x**9/9 + b**8*x**11/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_117():
    f = (a + b*x**2)**8/x**8
    F = -a**8/(7*x**7) - 8*a**7*b/(5*x**5) - 28*a**6*b**2/(3*x**3) - 56*a**5*b**3/x + 70*a**4*b**4*x + 56*a**3*b**5*x**3/3 + 28*a**2*b**6*x**5/5 + 8*a*b**7*x**7/7 + b**8*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_118():
    f = (a + b*x**2)**8/x**10
    F = -a**8/(9*x**9) - 8*a**7*b/(7*x**7) - 28*a**6*b**2/(5*x**5) - 56*a**5*b**3/(3*x**3) - 70*a**4*b**4/x + 56*a**3*b**5*x + 28*a**2*b**6*x**3/3 + 8*a*b**7*x**5/5 + b**8*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_119():
    f = (a + b*x**2)**8/x**12
    F = -a**8/(11*x**11) - 8*a**7*b/(9*x**9) - 4*a**6*b**2/x**7 - 56*a**5*b**3/(5*x**5) - 70*a**4*b**4/(3*x**3) - 56*a**3*b**5/x + 28*a**2*b**6*x + 8*a*b**7*x**3/3 + b**8*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_120():
    f = (a + b*x**2)**8/x**14
    F = -a**8/(13*x**13) - 8*a**7*b/(11*x**11) - 28*a**6*b**2/(9*x**9) - 8*a**5*b**3/x**7 - 14*a**4*b**4/x**5 - 56*a**3*b**5/(3*x**3) - 28*a**2*b**6/x + 8*a*b**7*x + b**8*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_121():
    f = (a + b*x**2)**8/x**16
    F = -a**8/(15*x**15) - 8*a**7*b/(13*x**13) - 28*a**6*b**2/(11*x**11) - 56*a**5*b**3/(9*x**9) - 10*a**4*b**4/x**7 - 56*a**3*b**5/(5*x**5) - 28*a**2*b**6/(3*x**3) - 8*a*b**7/x + b**8*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_122():
    f = (a + b*x**2)**8/x**18
    F = -a**8/(17*x**17) - 8*a**7*b/(15*x**15) - 28*a**6*b**2/(13*x**13) - 56*a**5*b**3/(11*x**11) - 70*a**4*b**4/(9*x**9) - 8*a**3*b**5/x**7 - 28*a**2*b**6/(5*x**5) - 8*a*b**7/(3*x**3) - b**8/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_123():
    f = (a + b*x**2)**8/x**20
    F = -a**8/(19*x**19) - 8*a**7*b/(17*x**17) - 28*a**6*b**2/(15*x**15) - 56*a**5*b**3/(13*x**13) - 70*a**4*b**4/(11*x**11) - 56*a**3*b**5/(9*x**9) - 4*a**2*b**6/x**7 - 8*a*b**7/(5*x**5) - b**8/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_124():
    f = x**11/(a + b*x**2)
    F = -a**5*log(a + b*x**2)/(2*b**6) + a**4*x**2/(2*b**5) - a**3*x**4/(4*b**4) + a**2*x**6/(6*b**3) - a*x**8/(8*b**2) + x**10/(10*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_125():
    f = x**10/(a + b*x**2)
    F = -a**(sympy.S(9)/2)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(11)/2) + a**4*x/b**5 - a**3*x**3/(3*b**4) + a**2*x**5/(5*b**3) - a*x**7/(7*b**2) + x**9/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_126():
    f = x**9/(a + b*x**2)
    F = a**4*log(a + b*x**2)/(2*b**5) - a**3*x**2/(2*b**4) + a**2*x**4/(4*b**3) - a*x**6/(6*b**2) + x**8/(8*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_127():
    f = x**8/(a + b*x**2)
    F = a**(sympy.S(7)/2)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(9)/2) - a**3*x/b**4 + a**2*x**3/(3*b**3) - a*x**5/(5*b**2) + x**7/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_128():
    f = x**7/(a + b*x**2)
    F = -a**3*log(a + b*x**2)/(2*b**4) + a**2*x**2/(2*b**3) - a*x**4/(4*b**2) + x**6/(6*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_129():
    f = x**6/(a + b*x**2)
    F = -a**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(7)/2) + a**2*x/b**3 - a*x**3/(3*b**2) + x**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_130():
    f = x**5/(a + b*x**2)
    F = a**2*log(a + b*x**2)/(2*b**3) - a*x**2/(2*b**2) + x**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_131():
    f = x**4/(a + b*x**2)
    F = a**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(5)/2) - a*x/b**2 + x**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_132():
    f = x**3/(a + b*x**2)
    F = -a*log(a + b*x**2)/(2*b**2) + x**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_133():
    f = x**2/(a + b*x**2)
    F = -sqrt(a)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(3)/2) + x/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_134():
    f = x/(a + b*x**2)
    F = log(a + b*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_135():
    f = 1/(a + b*x**2)
    F = atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_136():
    f = 1/(x*(a + b*x**2))
    F = log(x)/a - log(a + b*x**2)/(2*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_137():
    f = 1/(x**2*(a + b*x**2))
    F = -1/(a*x) - sqrt(b)*atan(sqrt(b)*x/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_138():
    f = 1/(x**3*(a + b*x**2))
    F = -1/(2*a*x**2) - b*log(x)/a**2 + b*log(a + b*x**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_139():
    f = 1/(x**4*(a + b*x**2))
    F = -1/(3*a*x**3) + b/(a**2*x) + b**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_140():
    f = 1/(x**5*(a + b*x**2))
    F = -1/(4*a*x**4) + b/(2*a**2*x**2) + b**2*log(x)/a**3 - b**2*log(a + b*x**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_141():
    f = 1/(x**6*(a + b*x**2))
    F = -1/(5*a*x**5) + b/(3*a**2*x**3) - b**2/(a**3*x) - b**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/a**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_142():
    f = 1/(x**7*(a + b*x**2))
    F = -1/(6*a*x**6) + b/(4*a**2*x**4) - b**2/(2*a**3*x**2) - b**3*log(x)/a**4 + b**3*log(a + b*x**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_143():
    f = 1/(x**8*(a + b*x**2))
    F = -1/(7*a*x**7) + b/(5*a**2*x**5) - b**2/(3*a**3*x**3) + b**3/(a**4*x) + b**(sympy.S(7)/2)*atan(sqrt(b)*x/sqrt(a))/a**(sympy.S(9)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_144():
    f = 1/(x**9*(a + b*x**2))
    F = -1/(8*a*x**8) + b/(6*a**2*x**6) - b**2/(4*a**3*x**4) + b**3/(2*a**4*x**2) + b**4*log(x)/a**5 - b**4*log(a + b*x**2)/(2*a**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_145():
    f = x**13/(a + b*x**2)**2
    F = -a**6/(2*b**7*(a + b*x**2)) - 3*a**5*log(a + b*x**2)/b**7 + 5*a**4*x**2/(2*b**6) - a**3*x**4/b**5 + a**2*x**6/(2*b**4) - a*x**8/(4*b**3) + x**10/(10*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_146():
    f = x**12/(a + b*x**2)**2
    F = -11*a**(sympy.S(9)/2)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(13)/2)) + 11*a**4*x/(2*b**6) - 11*a**3*x**3/(6*b**5) + 11*a**2*x**5/(10*b**4) - 11*a*x**7/(14*b**3) - x**11/(2*b*(a + b*x**2)) + 11*x**9/(18*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_147():
    f = x**11/(a + b*x**2)**2
    F = a**5/(2*b**6*(a + b*x**2)) + 5*a**4*log(a + b*x**2)/(2*b**6) - 2*a**3*x**2/b**5 + 3*a**2*x**4/(4*b**4) - a*x**6/(3*b**3) + x**8/(8*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_148():
    f = x**10/(a + b*x**2)**2
    F = 9*a**(sympy.S(7)/2)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(11)/2)) - 9*a**3*x/(2*b**5) + 3*a**2*x**3/(2*b**4) - 9*a*x**5/(10*b**3) - x**9/(2*b*(a + b*x**2)) + 9*x**7/(14*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_149():
    f = x**9/(a + b*x**2)**2
    F = -a**4/(2*b**5*(a + b*x**2)) - 2*a**3*log(a + b*x**2)/b**5 + 3*a**2*x**2/(2*b**4) - a*x**4/(2*b**3) + x**6/(6*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_150():
    f = x**8/(a + b*x**2)**2
    F = -7*a**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(9)/2)) + 7*a**2*x/(2*b**4) - 7*a*x**3/(6*b**3) - x**7/(2*b*(a + b*x**2)) + 7*x**5/(10*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_151():
    f = x**7/(a + b*x**2)**2
    F = a**3/(2*b**4*(a + b*x**2)) + 3*a**2*log(a + b*x**2)/(2*b**4) - a*x**2/b**3 + x**4/(4*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_152():
    f = x**6/(a + b*x**2)**2
    F = 5*a**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(7)/2)) - 5*a*x/(2*b**3) - x**5/(2*b*(a + b*x**2)) + 5*x**3/(6*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_153():
    f = x**5/(a + b*x**2)**2
    F = -a**2/(2*b**3*(a + b*x**2)) - a*log(a + b*x**2)/b**3 + x**2/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_154():
    f = x**4/(a + b*x**2)**2
    F = -3*sqrt(a)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(5)/2)) - x**3/(2*b*(a + b*x**2)) + 3*x/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_155():
    f = x**3/(a + b*x**2)**2
    F = a/(2*b**2*(a + b*x**2)) + log(a + b*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_156():
    f = x**2/(a + b*x**2)**2
    F = -x/(2*b*(a + b*x**2)) + atan(sqrt(b)*x/sqrt(a))/(2*sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_157():
    f = x/(a + b*x**2)**2
    F = -1/(2*b*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_158():
    f = (a + b*x**2)**(-2)
    F = x/(2*a*(a + b*x**2)) + atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_159():
    f = 1/(x*(a + b*x**2)**2)
    F = 1/(2*a*(a + b*x**2)) + log(x)/a**2 - log(a + b*x**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_160():
    f = 1/(x**2*(a + b*x**2)**2)
    F = 1/(2*a*x*(a + b*x**2)) - 3/(2*a**2*x) - 3*sqrt(b)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_161():
    f = 1/(x**3*(a + b*x**2)**2)
    F = -b/(2*a**2*(a + b*x**2)) - 1/(2*a**2*x**2) - 2*b*log(x)/a**3 + b*log(a + b*x**2)/a**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_162():
    f = 1/(x**4*(a + b*x**2)**2)
    F = 1/(2*a*x**3*(a + b*x**2)) - 5/(6*a**2*x**3) + 5*b/(2*a**3*x) + 5*b**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_163():
    f = 1/(x**5*(a + b*x**2)**2)
    F = -1/(4*a**2*x**4) + b**2/(2*a**3*(a + b*x**2)) + b/(a**3*x**2) + 3*b**2*log(x)/a**4 - 3*b**2*log(a + b*x**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_164():
    f = 1/(x**6*(a + b*x**2)**2)
    F = 1/(2*a*x**5*(a + b*x**2)) - 7/(10*a**2*x**5) + 7*b/(6*a**3*x**3) - 7*b**2/(2*a**4*x) - 7*b**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_165():
    f = 1/(x**7*(a + b*x**2)**2)
    F = -1/(6*a**2*x**6) + b/(2*a**3*x**4) - b**3/(2*a**4*(a + b*x**2)) - 3*b**2/(2*a**4*x**2) - 4*b**3*log(x)/a**5 + 2*b**3*log(a + b*x**2)/a**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_166():
    f = 1/(x**8*(a + b*x**2)**2)
    F = 1/(2*a*x**7*(a + b*x**2)) - 9/(14*a**2*x**7) + 9*b/(10*a**3*x**5) - 3*b**2/(2*a**4*x**3) + 9*b**3/(2*a**5*x) + 9*b**(sympy.S(7)/2)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_167():
    f = 1/(x**9*(a + b*x**2)**2)
    F = -1/(8*a**2*x**8) + b/(3*a**3*x**6) - 3*b**2/(4*a**4*x**4) + b**4/(2*a**5*(a + b*x**2)) + 2*b**3/(a**5*x**2) + 5*b**4*log(x)/a**6 - 5*b**4*log(a + b*x**2)/(2*a**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_168():
    f = x**15/(a + b*x**2)**3
    F = a**7/(4*b**8*(a + b*x**2)**2) - 7*a**6/(2*b**8*(a + b*x**2)) - 21*a**5*log(a + b*x**2)/(2*b**8) + 15*a**4*x**2/(2*b**7) - 5*a**3*x**4/(2*b**6) + a**2*x**6/b**5 - 3*a*x**8/(8*b**4) + x**10/(10*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_169():
    f = x**13/(a + b*x**2)**3
    F = -a**6/(4*b**7*(a + b*x**2)**2) + 3*a**5/(b**7*(a + b*x**2)) + 15*a**4*log(a + b*x**2)/(2*b**7) - 5*a**3*x**2/b**6 + 3*a**2*x**4/(2*b**5) - a*x**6/(2*b**4) + x**8/(8*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_170():
    f = x**11/(a + b*x**2)**3
    F = a**5/(4*b**6*(a + b*x**2)**2) - 5*a**4/(2*b**6*(a + b*x**2)) - 5*a**3*log(a + b*x**2)/b**6 + 3*a**2*x**2/b**5 - 3*a*x**4/(4*b**4) + x**6/(6*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_171():
    f = x**9/(a + b*x**2)**3
    F = -a**4/(4*b**5*(a + b*x**2)**2) + 2*a**3/(b**5*(a + b*x**2)) + 3*a**2*log(a + b*x**2)/b**5 - 3*a*x**2/(2*b**4) + x**4/(4*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_172():
    f = x**7/(a + b*x**2)**3
    F = a**3/(4*b**4*(a + b*x**2)**2) - 3*a**2/(2*b**4*(a + b*x**2)) - 3*a*log(a + b*x**2)/(2*b**4) + x**2/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_173():
    f = x**5/(a + b*x**2)**3
    F = -a**2/(4*b**3*(a + b*x**2)**2) + a/(b**3*(a + b*x**2)) + log(a + b*x**2)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_174():
    f = x**3/(a + b*x**2)**3
    F = x**4/(4*a*(a + b*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_175():
    f = x/(a + b*x**2)**3
    F = -1/(4*b*(a + b*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_176():
    f = 1/(x*(a + b*x**2)**3)
    F = 1/(4*a*(a + b*x**2)**2) + 1/(2*a**2*(a + b*x**2)) + log(x)/a**3 - log(a + b*x**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_177():
    f = 1/(x**3*(a + b*x**2)**3)
    F = -b/(4*a**2*(a + b*x**2)**2) - b/(a**3*(a + b*x**2)) - 1/(2*a**3*x**2) - 3*b*log(x)/a**4 + 3*b*log(a + b*x**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_178():
    f = 1/(x**5*(a + b*x**2)**3)
    F = b**2/(4*a**3*(a + b*x**2)**2) - 1/(4*a**3*x**4) + 3*b**2/(2*a**4*(a + b*x**2)) + 3*b/(2*a**4*x**2) + 6*b**2*log(x)/a**5 - 3*b**2*log(a + b*x**2)/a**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_179():
    f = 1/(x**7*(a + b*x**2)**3)
    F = -1/(6*a**3*x**6) - b**3/(4*a**4*(a + b*x**2)**2) + 3*b/(4*a**4*x**4) - 2*b**3/(a**5*(a + b*x**2)) - 3*b**2/(a**5*x**2) - 10*b**3*log(x)/a**6 + 5*b**3*log(a + b*x**2)/a**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_180():
    f = 1/(x**9*(a + b*x**2)**3)
    F = -1/(8*a**3*x**8) + b/(2*a**4*x**6) + b**4/(4*a**5*(a + b*x**2)**2) - 3*b**2/(2*a**5*x**4) + 5*b**4/(2*a**6*(a + b*x**2)) + 5*b**3/(a**6*x**2) + 15*b**4*log(x)/a**7 - 15*b**4*log(a + b*x**2)/(2*a**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_181():
    f = x**12/(a + b*x**2)**3
    F = 99*a**(sympy.S(7)/2)*atan(sqrt(b)*x/sqrt(a))/(8*b**(sympy.S(13)/2)) - 99*a**3*x/(8*b**6) + 33*a**2*x**3/(8*b**5) - 99*a*x**5/(40*b**4) - x**11/(4*b*(a + b*x**2)**2) - 11*x**9/(8*b**2*(a + b*x**2)) + 99*x**7/(56*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_182():
    f = x**10/(a + b*x**2)**3
    F = -63*a**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(8*b**(sympy.S(11)/2)) + 63*a**2*x/(8*b**5) - 21*a*x**3/(8*b**4) - x**9/(4*b*(a + b*x**2)**2) - 9*x**7/(8*b**2*(a + b*x**2)) + 63*x**5/(40*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_183():
    f = x**8/(a + b*x**2)**3
    F = 35*a**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(8*b**(sympy.S(9)/2)) - 35*a*x/(8*b**4) - x**7/(4*b*(a + b*x**2)**2) - 7*x**5/(8*b**2*(a + b*x**2)) + 35*x**3/(24*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_184():
    f = x**6/(a + b*x**2)**3
    F = -15*sqrt(a)*atan(sqrt(b)*x/sqrt(a))/(8*b**(sympy.S(7)/2)) - x**5/(4*b*(a + b*x**2)**2) - 5*x**3/(8*b**2*(a + b*x**2)) + 15*x/(8*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_185():
    f = x**4/(a + b*x**2)**3
    F = -x**3/(4*b*(a + b*x**2)**2) - 3*x/(8*b**2*(a + b*x**2)) + 3*atan(sqrt(b)*x/sqrt(a))/(8*sqrt(a)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_186():
    f = x**2/(a + b*x**2)**3
    F = -x/(4*b*(a + b*x**2)**2) + x/(8*a*b*(a + b*x**2)) + atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(3)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_187():
    f = (a + b*x**2)**(-3)
    F = x/(4*a*(a + b*x**2)**2) + 3*x/(8*a**2*(a + b*x**2)) + 3*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_188():
    f = 1/(x**2*(a + b*x**2)**3)
    F = 1/(4*a*x*(a + b*x**2)**2) + 5/(8*a**2*x*(a + b*x**2)) - 15/(8*a**3*x) - 15*sqrt(b)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_189():
    f = 1/(x**4*(a + b*x**2)**3)
    F = 1/(4*a*x**3*(a + b*x**2)**2) + 7/(8*a**2*x**3*(a + b*x**2)) - 35/(24*a**3*x**3) + 35*b/(8*a**4*x) + 35*b**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_190():
    f = 1/(x**6*(a + b*x**2)**3)
    F = 1/(4*a*x**5*(a + b*x**2)**2) + 9/(8*a**2*x**5*(a + b*x**2)) - 63/(40*a**3*x**5) + 21*b/(8*a**4*x**3) - 63*b**2/(8*a**5*x) - 63*b**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_191():
    f = 1/(x**8*(a + b*x**2)**3)
    F = 1/(4*a*x**7*(a + b*x**2)**2) + 11/(8*a**2*x**7*(a + b*x**2)) - 99/(56*a**3*x**7) + 99*b/(40*a**4*x**5) - 33*b**2/(8*a**5*x**3) + 99*b**3/(8*a**6*x) + 99*b**(sympy.S(7)/2)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_192():
    f = x**25/(a + b*x**2)**10
    F = -a**12/(18*b**13*(a + b*x**2)**9) + 3*a**11/(4*b**13*(a + b*x**2)**8) - 33*a**10/(7*b**13*(a + b*x**2)**7) + 55*a**9/(3*b**13*(a + b*x**2)**6) - 99*a**8/(2*b**13*(a + b*x**2)**5) + 99*a**7/(b**13*(a + b*x**2)**4) - 154*a**6/(b**13*(a + b*x**2)**3) + 198*a**5/(b**13*(a + b*x**2)**2) - 495*a**4/(2*b**13*(a + b*x**2)) - 110*a**3*log(a + b*x**2)/b**13 + 55*a**2*x**2/(2*b**12) - 5*a*x**4/(2*b**11) + x**6/(6*b**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_193():
    f = x**23/(a + b*x**2)**10
    F = a**11/(18*b**12*(a + b*x**2)**9) - 11*a**10/(16*b**12*(a + b*x**2)**8) + 55*a**9/(14*b**12*(a + b*x**2)**7) - 55*a**8/(4*b**12*(a + b*x**2)**6) + 33*a**7/(b**12*(a + b*x**2)**5) - 231*a**6/(4*b**12*(a + b*x**2)**4) + 77*a**5/(b**12*(a + b*x**2)**3) - 165*a**4/(2*b**12*(a + b*x**2)**2) + 165*a**3/(2*b**12*(a + b*x**2)) + 55*a**2*log(a + b*x**2)/(2*b**12) - 5*a*x**2/b**11 + x**4/(4*b**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_194():
    f = x**21/(a + b*x**2)**10
    F = -a**10/(18*b**11*(a + b*x**2)**9) + 5*a**9/(8*b**11*(a + b*x**2)**8) - 45*a**8/(14*b**11*(a + b*x**2)**7) + 10*a**7/(b**11*(a + b*x**2)**6) - 21*a**6/(b**11*(a + b*x**2)**5) + 63*a**5/(2*b**11*(a + b*x**2)**4) - 35*a**4/(b**11*(a + b*x**2)**3) + 30*a**3/(b**11*(a + b*x**2)**2) - 45*a**2/(2*b**11*(a + b*x**2)) - 5*a*log(a + b*x**2)/b**11 + x**2/(2*b**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_195():
    f = x**19/(a + b*x**2)**10
    F = a**9/(18*b**10*(a + b*x**2)**9) - 9*a**8/(16*b**10*(a + b*x**2)**8) + 18*a**7/(7*b**10*(a + b*x**2)**7) - 7*a**6/(b**10*(a + b*x**2)**6) + 63*a**5/(5*b**10*(a + b*x**2)**5) - 63*a**4/(4*b**10*(a + b*x**2)**4) + 14*a**3/(b**10*(a + b*x**2)**3) - 9*a**2/(b**10*(a + b*x**2)**2) + 9*a/(2*b**10*(a + b*x**2)) + log(a + b*x**2)/(2*b**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_196():
    f = x**17/(a + b*x**2)**10
    F = x**18/(18*a*(a + b*x**2)**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_197():
    f = x**15/(a + b*x**2)**10
    F = x**16/(18*a*(a + b*x**2)**9) + x**16/(144*a**2*(a + b*x**2)**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_198():
    f = x**13/(a + b*x**2)**10
    F = x**14/(18*a*(a + b*x**2)**9) + x**14/(72*a**2*(a + b*x**2)**8) + x**14/(504*a**3*(a + b*x**2)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_199():
    f = x**11/(a + b*x**2)**10
    F = x**12/(18*a*(a + b*x**2)**9) + x**12/(48*a**2*(a + b*x**2)**8) + x**12/(168*a**3*(a + b*x**2)**7) + x**12/(1008*a**4*(a + b*x**2)**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_200():
    f = x**9/(a + b*x**2)**10
    F = -a**4/(18*b**5*(a + b*x**2)**9) + a**3/(4*b**5*(a + b*x**2)**8) - 3*a**2/(7*b**5*(a + b*x**2)**7) + a/(3*b**5*(a + b*x**2)**6) - 1/(10*b**5*(a + b*x**2)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_201():
    f = x**7/(a + b*x**2)**10
    F = a**3/(18*b**4*(a + b*x**2)**9) - 3*a**2/(16*b**4*(a + b*x**2)**8) + 3*a/(14*b**4*(a + b*x**2)**7) - 1/(12*b**4*(a + b*x**2)**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_202():
    f = x**5/(a + b*x**2)**10
    F = -a**2/(18*b**3*(a + b*x**2)**9) + a/(8*b**3*(a + b*x**2)**8) - 1/(14*b**3*(a + b*x**2)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_203():
    f = x**3/(a + b*x**2)**10
    F = a/(18*b**2*(a + b*x**2)**9) - 1/(16*b**2*(a + b*x**2)**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_204():
    f = x/(a + b*x**2)**10
    F = -1/(18*b*(a + b*x**2)**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_205():
    f = 1/(x*(a + b*x**2)**10)
    F = 1/(18*a*(a + b*x**2)**9) + 1/(16*a**2*(a + b*x**2)**8) + 1/(14*a**3*(a + b*x**2)**7) + 1/(12*a**4*(a + b*x**2)**6) + 1/(10*a**5*(a + b*x**2)**5) + 1/(8*a**6*(a + b*x**2)**4) + 1/(6*a**7*(a + b*x**2)**3) + 1/(4*a**8*(a + b*x**2)**2) + 1/(2*a**9*(a + b*x**2)) + log(x)/a**10 - log(a + b*x**2)/(2*a**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_206():
    f = 1/(x**3*(a + b*x**2)**10)
    F = -b/(18*a**2*(a + b*x**2)**9) - b/(8*a**3*(a + b*x**2)**8) - 3*b/(14*a**4*(a + b*x**2)**7) - b/(3*a**5*(a + b*x**2)**6) - b/(2*a**6*(a + b*x**2)**5) - 3*b/(4*a**7*(a + b*x**2)**4) - 7*b/(6*a**8*(a + b*x**2)**3) - 2*b/(a**9*(a + b*x**2)**2) - 9*b/(2*a**10*(a + b*x**2)) - 1/(2*a**10*x**2) - 10*b*log(x)/a**11 + 5*b*log(a + b*x**2)/a**11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_207():
    f = 1/(x**5*(a + b*x**2)**10)
    F = b**2/(18*a**3*(a + b*x**2)**9) + 3*b**2/(16*a**4*(a + b*x**2)**8) + 3*b**2/(7*a**5*(a + b*x**2)**7) + 5*b**2/(6*a**6*(a + b*x**2)**6) + 3*b**2/(2*a**7*(a + b*x**2)**5) + 21*b**2/(8*a**8*(a + b*x**2)**4) + 14*b**2/(3*a**9*(a + b*x**2)**3) + 9*b**2/(a**10*(a + b*x**2)**2) - 1/(4*a**10*x**4) + 45*b**2/(2*a**11*(a + b*x**2)) + 5*b/(a**11*x**2) + 55*b**2*log(x)/a**12 - 55*b**2*log(a + b*x**2)/(2*a**12)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_208():
    f = 1/(x**7*(a + b*x**2)**10)
    F = -b**3/(18*a**4*(a + b*x**2)**9) - b**3/(4*a**5*(a + b*x**2)**8) - 5*b**3/(7*a**6*(a + b*x**2)**7) - 5*b**3/(3*a**7*(a + b*x**2)**6) - 7*b**3/(2*a**8*(a + b*x**2)**5) - 7*b**3/(a**9*(a + b*x**2)**4) - 14*b**3/(a**10*(a + b*x**2)**3) - 1/(6*a**10*x**6) - 30*b**3/(a**11*(a + b*x**2)**2) + 5*b/(2*a**11*x**4) - 165*b**3/(2*a**12*(a + b*x**2)) - 55*b**2/(2*a**12*x**2) - 220*b**3*log(x)/a**13 + 110*b**3*log(a + b*x**2)/a**13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_209():
    f = x**24/(a + b*x**2)**10
    F = -7436429*a**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(65536*b**(sympy.S(25)/2)) + 7436429*a**2*x/(65536*b**12) - 7436429*a*x**3/(196608*b**11) - x**23/(18*b*(a + b*x**2)**9) - 23*x**21/(288*b**2*(a + b*x**2)**8) - 23*x**19/(192*b**3*(a + b*x**2)**7) - 437*x**17/(2304*b**4*(a + b*x**2)**6) - 7429*x**15/(23040*b**5*(a + b*x**2)**5) - 7429*x**13/(12288*b**6*(a + b*x**2)**4) - 96577*x**11/(73728*b**7*(a + b*x**2)**3) - 1062347*x**9/(294912*b**8*(a + b*x**2)**2) - 1062347*x**7/(65536*b**9*(a + b*x**2)) + 7436429*x**5/(327680*b**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_210():
    f = x**22/(a + b*x**2)**10
    F = 1616615*a**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(65536*b**(sympy.S(23)/2)) - 1616615*a*x/(65536*b**11) - x**21/(18*b*(a + b*x**2)**9) - 7*x**19/(96*b**2*(a + b*x**2)**8) - 19*x**17/(192*b**3*(a + b*x**2)**7) - 323*x**15/(2304*b**4*(a + b*x**2)**6) - 323*x**13/(1536*b**5*(a + b*x**2)**5) - 4199*x**11/(12288*b**6*(a + b*x**2)**4) - 46189*x**9/(73728*b**7*(a + b*x**2)**3) - 46189*x**7/(32768*b**8*(a + b*x**2)**2) - 323323*x**5/(65536*b**9*(a + b*x**2)) + 1616615*x**3/(196608*b**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_211():
    f = x**20/(a + b*x**2)**10
    F = -230945*sqrt(a)*atan(sqrt(b)*x/sqrt(a))/(65536*b**(sympy.S(21)/2)) - x**19/(18*b*(a + b*x**2)**9) - 19*x**17/(288*b**2*(a + b*x**2)**8) - 323*x**15/(4032*b**3*(a + b*x**2)**7) - 1615*x**13/(16128*b**4*(a + b*x**2)**6) - 4199*x**11/(32256*b**5*(a + b*x**2)**5) - 46189*x**9/(258048*b**6*(a + b*x**2)**4) - 46189*x**7/(172032*b**7*(a + b*x**2)**3) - 46189*x**5/(98304*b**8*(a + b*x**2)**2) - 230945*x**3/(196608*b**9*(a + b*x**2)) + 230945*x/(65536*b**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_212():
    f = x**18/(a + b*x**2)**10
    F = -x**17/(18*b*(a + b*x**2)**9) - 17*x**15/(288*b**2*(a + b*x**2)**8) - 85*x**13/(1344*b**3*(a + b*x**2)**7) - 1105*x**11/(16128*b**4*(a + b*x**2)**6) - 2431*x**9/(32256*b**5*(a + b*x**2)**5) - 2431*x**7/(28672*b**6*(a + b*x**2)**4) - 2431*x**5/(24576*b**7*(a + b*x**2)**3) - 12155*x**3/(98304*b**8*(a + b*x**2)**2) - 12155*x/(65536*b**9*(a + b*x**2)) + 12155*atan(sqrt(b)*x/sqrt(a))/(65536*sqrt(a)*b**(sympy.S(19)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_213():
    f = x**16/(a + b*x**2)**10
    F = -x**15/(18*b*(a + b*x**2)**9) - 5*x**13/(96*b**2*(a + b*x**2)**8) - 65*x**11/(1344*b**3*(a + b*x**2)**7) - 715*x**9/(16128*b**4*(a + b*x**2)**6) - 143*x**7/(3584*b**5*(a + b*x**2)**5) - 143*x**5/(4096*b**6*(a + b*x**2)**4) - 715*x**3/(24576*b**7*(a + b*x**2)**3) - 715*x/(32768*b**8*(a + b*x**2)**2) + 715*x/(65536*a*b**8*(a + b*x**2)) + 715*atan(sqrt(b)*x/sqrt(a))/(65536*a**(sympy.S(3)/2)*b**(sympy.S(17)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_214():
    f = x**14/(a + b*x**2)**10
    F = -x**13/(18*b*(a + b*x**2)**9) - 13*x**11/(288*b**2*(a + b*x**2)**8) - 143*x**9/(4032*b**3*(a + b*x**2)**7) - 143*x**7/(5376*b**4*(a + b*x**2)**6) - 143*x**5/(7680*b**5*(a + b*x**2)**5) - 143*x**3/(12288*b**6*(a + b*x**2)**4) - 143*x/(24576*b**7*(a + b*x**2)**3) + 143*x/(98304*a*b**7*(a + b*x**2)**2) + 143*x/(65536*a**2*b**7*(a + b*x**2)) + 143*atan(sqrt(b)*x/sqrt(a))/(65536*a**(sympy.S(5)/2)*b**(sympy.S(15)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_215():
    f = x**12/(a + b*x**2)**10
    F = -x**11/(18*b*(a + b*x**2)**9) - 11*x**9/(288*b**2*(a + b*x**2)**8) - 11*x**7/(448*b**3*(a + b*x**2)**7) - 11*x**5/(768*b**4*(a + b*x**2)**6) - 11*x**3/(1536*b**5*(a + b*x**2)**5) - 11*x/(4096*b**6*(a + b*x**2)**4) + 11*x/(24576*a*b**6*(a + b*x**2)**3) + 55*x/(98304*a**2*b**6*(a + b*x**2)**2) + 55*x/(65536*a**3*b**6*(a + b*x**2)) + 55*atan(sqrt(b)*x/sqrt(a))/(65536*a**(sympy.S(7)/2)*b**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_216():
    f = x**10/(a + b*x**2)**10
    F = -x**9/(18*b*(a + b*x**2)**9) - x**7/(32*b**2*(a + b*x**2)**8) - x**5/(64*b**3*(a + b*x**2)**7) - 5*x**3/(768*b**4*(a + b*x**2)**6) - x/(512*b**5*(a + b*x**2)**5) + x/(4096*a*b**5*(a + b*x**2)**4) + 7*x/(24576*a**2*b**5*(a + b*x**2)**3) + 35*x/(98304*a**3*b**5*(a + b*x**2)**2) + 35*x/(65536*a**4*b**5*(a + b*x**2)) + 35*atan(sqrt(b)*x/sqrt(a))/(65536*a**(sympy.S(9)/2)*b**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_217():
    f = x**8/(a + b*x**2)**10
    F = -x**7/(18*b*(a + b*x**2)**9) - 7*x**5/(288*b**2*(a + b*x**2)**8) - 5*x**3/(576*b**3*(a + b*x**2)**7) - 5*x/(2304*b**4*(a + b*x**2)**6) + x/(4608*a*b**4*(a + b*x**2)**5) + x/(4096*a**2*b**4*(a + b*x**2)**4) + 7*x/(24576*a**3*b**4*(a + b*x**2)**3) + 35*x/(98304*a**4*b**4*(a + b*x**2)**2) + 35*x/(65536*a**5*b**4*(a + b*x**2)) + 35*atan(sqrt(b)*x/sqrt(a))/(65536*a**(sympy.S(11)/2)*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_218():
    f = x**6/(a + b*x**2)**10
    F = -x**5/(18*b*(a + b*x**2)**9) - 5*x**3/(288*b**2*(a + b*x**2)**8) - 5*x/(1344*b**3*(a + b*x**2)**7) + 5*x/(16128*a*b**3*(a + b*x**2)**6) + 11*x/(32256*a**2*b**3*(a + b*x**2)**5) + 11*x/(28672*a**3*b**3*(a + b*x**2)**4) + 11*x/(24576*a**4*b**3*(a + b*x**2)**3) + 55*x/(98304*a**5*b**3*(a + b*x**2)**2) + 55*x/(65536*a**6*b**3*(a + b*x**2)) + 55*atan(sqrt(b)*x/sqrt(a))/(65536*a**(sympy.S(13)/2)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_219():
    f = x**4/(a + b*x**2)**10
    F = -x**3/(18*b*(a + b*x**2)**9) - x/(96*b**2*(a + b*x**2)**8) + x/(1344*a*b**2*(a + b*x**2)**7) + 13*x/(16128*a**2*b**2*(a + b*x**2)**6) + 143*x/(161280*a**3*b**2*(a + b*x**2)**5) + 143*x/(143360*a**4*b**2*(a + b*x**2)**4) + 143*x/(122880*a**5*b**2*(a + b*x**2)**3) + 143*x/(98304*a**6*b**2*(a + b*x**2)**2) + 143*x/(65536*a**7*b**2*(a + b*x**2)) + 143*atan(sqrt(b)*x/sqrt(a))/(65536*a**(sympy.S(15)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_220():
    f = x**2/(a + b*x**2)**10
    F = -x/(18*b*(a + b*x**2)**9) + x/(288*a*b*(a + b*x**2)**8) + 5*x/(1344*a**2*b*(a + b*x**2)**7) + 65*x/(16128*a**3*b*(a + b*x**2)**6) + 143*x/(32256*a**4*b*(a + b*x**2)**5) + 143*x/(28672*a**5*b*(a + b*x**2)**4) + 143*x/(24576*a**6*b*(a + b*x**2)**3) + 715*x/(98304*a**7*b*(a + b*x**2)**2) + 715*x/(65536*a**8*b*(a + b*x**2)) + 715*atan(sqrt(b)*x/sqrt(a))/(65536*a**(sympy.S(17)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_221():
    f = (a + b*x**2)**(-10)
    F = x/(18*a*(a + b*x**2)**9) + 17*x/(288*a**2*(a + b*x**2)**8) + 85*x/(1344*a**3*(a + b*x**2)**7) + 1105*x/(16128*a**4*(a + b*x**2)**6) + 2431*x/(32256*a**5*(a + b*x**2)**5) + 2431*x/(28672*a**6*(a + b*x**2)**4) + 2431*x/(24576*a**7*(a + b*x**2)**3) + 12155*x/(98304*a**8*(a + b*x**2)**2) + 12155*x/(65536*a**9*(a + b*x**2)) + 12155*atan(sqrt(b)*x/sqrt(a))/(65536*a**(sympy.S(19)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_222():
    f = 1/(x**2*(a + b*x**2)**10)
    F = 1/(18*a*x*(a + b*x**2)**9) + 19/(288*a**2*x*(a + b*x**2)**8) + 323/(4032*a**3*x*(a + b*x**2)**7) + 1615/(16128*a**4*x*(a + b*x**2)**6) + 4199/(32256*a**5*x*(a + b*x**2)**5) + 46189/(258048*a**6*x*(a + b*x**2)**4) + 46189/(172032*a**7*x*(a + b*x**2)**3) + 46189/(98304*a**8*x*(a + b*x**2)**2) + 230945/(196608*a**9*x*(a + b*x**2)) - 230945/(65536*a**10*x) - 230945*sqrt(b)*atan(sqrt(b)*x/sqrt(a))/(65536*a**(sympy.S(21)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_223():
    f = 1/(x**4*(a + b*x**2)**10)
    F = 1/(18*a*x**3*(a + b*x**2)**9) + 7/(96*a**2*x**3*(a + b*x**2)**8) + 19/(192*a**3*x**3*(a + b*x**2)**7) + 323/(2304*a**4*x**3*(a + b*x**2)**6) + 323/(1536*a**5*x**3*(a + b*x**2)**5) + 4199/(12288*a**6*x**3*(a + b*x**2)**4) + 46189/(73728*a**7*x**3*(a + b*x**2)**3) + 46189/(32768*a**8*x**3*(a + b*x**2)**2) + 323323/(65536*a**9*x**3*(a + b*x**2)) - 1616615/(196608*a**10*x**3) + 1616615*b/(65536*a**11*x) + 1616615*b**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(65536*a**(sympy.S(23)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_224():
    f = 1/(x**6*(a + b*x**2)**10)
    F = 1/(18*a*x**5*(a + b*x**2)**9) + 23/(288*a**2*x**5*(a + b*x**2)**8) + 23/(192*a**3*x**5*(a + b*x**2)**7) + 437/(2304*a**4*x**5*(a + b*x**2)**6) + 7429/(23040*a**5*x**5*(a + b*x**2)**5) + 7429/(12288*a**6*x**5*(a + b*x**2)**4) + 96577/(73728*a**7*x**5*(a + b*x**2)**3) + 1062347/(294912*a**8*x**5*(a + b*x**2)**2) + 1062347/(65536*a**9*x**5*(a + b*x**2)) - 7436429/(327680*a**10*x**5) + 7436429*b/(196608*a**11*x**3) - 7436429*b**2/(65536*a**12*x) - 7436429*b**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(65536*a**(sympy.S(25)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_225():
    f = x**3/(a - b*x**2)
    F = -a*log(a - b*x**2)/(2*b**2) - x**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_226():
    f = x**2/(a - b*x**2)
    F = sqrt(a)*atanh(sqrt(b)*x/sqrt(a))/b**(sympy.S(3)/2) - x/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_227():
    f = x/(a - b*x**2)
    F = -log(a - b*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_228():
    f = 1/(a - b*x**2)
    F = atanh(sqrt(b)*x/sqrt(a))/(sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_229():
    f = 1/(x*(a - b*x**2))
    F = log(x)/a - log(a - b*x**2)/(2*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_230():
    f = 1/(x**2*(a - b*x**2))
    F = -1/(a*x) + sqrt(b)*atanh(sqrt(b)*x/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_231():
    f = 1/(x**3*(a - b*x**2))
    F = -1/(2*a*x**2) + b*log(x)/a**2 - b*log(a - b*x**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_232():
    f = x**3/(a - b*x**2)**2
    F = a/(2*b**2*(a - b*x**2)) + log(a - b*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_233():
    f = x**2/(a - b*x**2)**2
    F = x/(2*b*(a - b*x**2)) - atanh(sqrt(b)*x/sqrt(a))/(2*sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_234():
    f = x/(a - b*x**2)**2
    F = 1/(2*b*(a - b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_235():
    f = (a - b*x**2)**(-2)
    F = x/(2*a*(a - b*x**2)) + atanh(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_236():
    f = 1/(x*(a - b*x**2)**2)
    F = 1/(2*a*(a - b*x**2)) + log(x)/a**2 - log(a - b*x**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_237():
    f = 1/(x**2*(a - b*x**2)**2)
    F = 1/(2*a*x*(a - b*x**2)) - 3/(2*a**2*x) + 3*sqrt(b)*atanh(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_238():
    f = 1/(x**3*(a - b*x**2)**2)
    F = b/(2*a**2*(a - b*x**2)) - 1/(2*a**2*x**2) + 2*b*log(x)/a**3 - b*log(a - b*x**2)/a**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_239():
    f = x**3/(a - b*x**2)**3
    F = x**4/(4*a*(a - b*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_240():
    f = x**2/(a - b*x**2)**3
    F = x/(4*b*(a - b*x**2)**2) - x/(8*a*b*(a - b*x**2)) - atanh(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(3)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_241():
    f = x/(a - b*x**2)**3
    F = 1/(4*b*(a - b*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_242():
    f = (a - b*x**2)**(-3)
    F = x/(4*a*(a - b*x**2)**2) + 3*x/(8*a**2*(a - b*x**2)) + 3*atanh(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_243():
    f = 1/(x*(a - b*x**2)**3)
    F = 1/(4*a*(a - b*x**2)**2) + 1/(2*a**2*(a - b*x**2)) + log(x)/a**3 - log(a - b*x**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_244():
    f = 1/(x**2*(a - b*x**2)**3)
    F = 1/(4*a*x*(a - b*x**2)**2) + 5/(8*a**2*x*(a - b*x**2)) - 15/(8*a**3*x) + 15*sqrt(b)*atanh(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_245():
    f = 1/(x**3*(a - b*x**2)**3)
    F = b/(4*a**2*(a - b*x**2)**2) + b/(a**3*(a - b*x**2)) - 1/(2*a**3*x**2) + 3*b*log(x)/a**4 - 3*b*log(a - b*x**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_246():
    f = x**3/(a - b*x**2)**5
    F = a/(8*b**2*(a - b*x**2)**4) - 1/(6*b**2*(a - b*x**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_247():
    f = x**2/(a - b*x**2)**5
    F = x/(8*b*(a - b*x**2)**4) - x/(48*a*b*(a - b*x**2)**3) - 5*x/(192*a**2*b*(a - b*x**2)**2) - 5*x/(128*a**3*b*(a - b*x**2)) - 5*atanh(sqrt(b)*x/sqrt(a))/(128*a**(sympy.S(7)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_248():
    f = x/(a - b*x**2)**5
    F = 1/(8*b*(a - b*x**2)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_249():
    f = (a - b*x**2)**(-5)
    F = x/(8*a*(a - b*x**2)**4) + 7*x/(48*a**2*(a - b*x**2)**3) + 35*x/(192*a**3*(a - b*x**2)**2) + 35*x/(128*a**4*(a - b*x**2)) + 35*atanh(sqrt(b)*x/sqrt(a))/(128*a**(sympy.S(9)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_250():
    f = 1/(x*(a - b*x**2)**5)
    F = 1/(8*a*(a - b*x**2)**4) + 1/(6*a**2*(a - b*x**2)**3) + 1/(4*a**3*(a - b*x**2)**2) + 1/(2*a**4*(a - b*x**2)) + log(x)/a**5 - log(a - b*x**2)/(2*a**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_251():
    f = 1/(x**2*(a - b*x**2)**5)
    F = 1/(8*a*x*(a - b*x**2)**4) + 3/(16*a**2*x*(a - b*x**2)**3) + 21/(64*a**3*x*(a - b*x**2)**2) + 105/(128*a**4*x*(a - b*x**2)) - 315/(128*a**5*x) + 315*sqrt(b)*atanh(sqrt(b)*x/sqrt(a))/(128*a**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_252():
    f = 1/(x**3*(a - b*x**2)**5)
    F = b/(8*a**2*(a - b*x**2)**4) + b/(3*a**3*(a - b*x**2)**3) + 3*b/(4*a**4*(a - b*x**2)**2) + 2*b/(a**5*(a - b*x**2)) - 1/(2*a**5*x**2) + 5*b*log(x)/a**6 - 5*b*log(a - b*x**2)/(2*a**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_253():
    f = 1/(x*(b*x**2 + 1))
    F = log(x) - log(b*x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_254():
    f = 1/(x*(b*x**2 - 1))
    F = -log(x) + log(-b*x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_255():
    f = 1/(x**3*(b*x**2 + 1))
    F = -b*log(x) + b*log(b*x**2 + 1)/2 - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_256():
    f = 1/(x**3*(b*x**2 - 1))
    F = -b*log(x) + b*log(-b*x**2 + 1)/2 + 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_257():
    f = 1/(a*x**2 + a - 1)
    F = -atanh(sqrt(a)*x/sqrt(1 - a))/sqrt(a*(1 - a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_258():
    f = 1/(-c - d + x**2*(c - d))
    F = -atanh(x*sqrt(c - d)/sqrt(c + d))/(sqrt(c - d)*sqrt(c + d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_259():
    f = 1/(x*(b*x**2 + 1)**2)
    F = log(x) - log(b*x**2 + 1)/2 + 1/(2*b*x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_260():
    f = 1/(x*(b*x**2 - 1)**2)
    F = log(x) - log(-b*x**2 + 1)/2 + 1/(-2*b*x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_261():
    f = 1/(a + x**2*(-a*c + b))
    F = atan(x*sqrt(-a*c + b)/sqrt(a))/(sqrt(a)*sqrt(-a*c + b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_262():
    f = 1/(a - x**2*(-a*c + b))
    F = atanh(x*sqrt(-a*c + b)/sqrt(a))/(sqrt(a)*sqrt(-a*c + b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_263():
    f = 1/(c*(a - d) - x**2*(b - c))
    F = atanh(x*sqrt(b - c)/(sqrt(c)*sqrt(a - d)))/(sqrt(c)*sqrt(a - d)*sqrt(b - c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_264():
    f = x**(sympy.S(7)/2)*(a + b*x**2)
    F = 2*a*x**(sympy.S(9)/2)/9 + 2*b*x**(sympy.S(13)/2)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_265():
    f = x**(sympy.S(5)/2)*(a + b*x**2)
    F = 2*a*x**(sympy.S(7)/2)/7 + 2*b*x**(sympy.S(11)/2)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_266():
    f = x**(sympy.S(3)/2)*(a + b*x**2)
    F = 2*a*x**(sympy.S(5)/2)/5 + 2*b*x**(sympy.S(9)/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_267():
    f = sqrt(x)*(a + b*x**2)
    F = 2*a*x**(sympy.S(3)/2)/3 + 2*b*x**(sympy.S(7)/2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_268():
    f = (a + b*x**2)/sqrt(x)
    F = 2*a*sqrt(x) + 2*b*x**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_269():
    f = (a + b*x**2)/x**(sympy.S(3)/2)
    F = -2*a/sqrt(x) + 2*b*x**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_270():
    f = (a + b*x**2)/x**(sympy.S(5)/2)
    F = -2*a/(3*x**(sympy.S(3)/2)) + 2*b*sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_271():
    f = (a + b*x**2)/x**(sympy.S(7)/2)
    F = -2*a/(5*x**(sympy.S(5)/2)) - 2*b/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_272():
    f = x**(sympy.S(7)/2)*(a + b*x**2)**2
    F = 2*a**2*x**(sympy.S(9)/2)/9 + 4*a*b*x**(sympy.S(13)/2)/13 + 2*b**2*x**(sympy.S(17)/2)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_273():
    f = x**(sympy.S(5)/2)*(a + b*x**2)**2
    F = 2*a**2*x**(sympy.S(7)/2)/7 + 4*a*b*x**(sympy.S(11)/2)/11 + 2*b**2*x**(sympy.S(15)/2)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_274():
    f = x**(sympy.S(3)/2)*(a + b*x**2)**2
    F = 2*a**2*x**(sympy.S(5)/2)/5 + 4*a*b*x**(sympy.S(9)/2)/9 + 2*b**2*x**(sympy.S(13)/2)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_275():
    f = sqrt(x)*(a + b*x**2)**2
    F = 2*a**2*x**(sympy.S(3)/2)/3 + 4*a*b*x**(sympy.S(7)/2)/7 + 2*b**2*x**(sympy.S(11)/2)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_276():
    f = (a + b*x**2)**2/sqrt(x)
    F = 2*a**2*sqrt(x) + 4*a*b*x**(sympy.S(5)/2)/5 + 2*b**2*x**(sympy.S(9)/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_277():
    f = (a + b*x**2)**2/x**(sympy.S(3)/2)
    F = -2*a**2/sqrt(x) + 4*a*b*x**(sympy.S(3)/2)/3 + 2*b**2*x**(sympy.S(7)/2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_278():
    f = (a + b*x**2)**2/x**(sympy.S(5)/2)
    F = -2*a**2/(3*x**(sympy.S(3)/2)) + 4*a*b*sqrt(x) + 2*b**2*x**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_279():
    f = (a + b*x**2)**2/x**(sympy.S(7)/2)
    F = -2*a**2/(5*x**(sympy.S(5)/2)) - 4*a*b/sqrt(x) + 2*b**2*x**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_280():
    f = x**(sympy.S(7)/2)*(a + b*x**2)**3
    F = 2*a**3*x**(sympy.S(9)/2)/9 + 6*a**2*b*x**(sympy.S(13)/2)/13 + 6*a*b**2*x**(sympy.S(17)/2)/17 + 2*b**3*x**(sympy.S(21)/2)/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_281():
    f = x**(sympy.S(5)/2)*(a + b*x**2)**3
    F = 2*a**3*x**(sympy.S(7)/2)/7 + 6*a**2*b*x**(sympy.S(11)/2)/11 + 2*a*b**2*x**(sympy.S(15)/2)/5 + 2*b**3*x**(sympy.S(19)/2)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_282():
    f = x**(sympy.S(3)/2)*(a + b*x**2)**3
    F = 2*a**3*x**(sympy.S(5)/2)/5 + 2*a**2*b*x**(sympy.S(9)/2)/3 + 6*a*b**2*x**(sympy.S(13)/2)/13 + 2*b**3*x**(sympy.S(17)/2)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_283():
    f = sqrt(x)*(a + b*x**2)**3
    F = 2*a**3*x**(sympy.S(3)/2)/3 + 6*a**2*b*x**(sympy.S(7)/2)/7 + 6*a*b**2*x**(sympy.S(11)/2)/11 + 2*b**3*x**(sympy.S(15)/2)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_284():
    f = (a + b*x**2)**3/sqrt(x)
    F = 2*a**3*sqrt(x) + 6*a**2*b*x**(sympy.S(5)/2)/5 + 2*a*b**2*x**(sympy.S(9)/2)/3 + 2*b**3*x**(sympy.S(13)/2)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_285():
    f = (a + b*x**2)**3/x**(sympy.S(3)/2)
    F = -2*a**3/sqrt(x) + 2*a**2*b*x**(sympy.S(3)/2) + 6*a*b**2*x**(sympy.S(7)/2)/7 + 2*b**3*x**(sympy.S(11)/2)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_286():
    f = (a + b*x**2)**3/x**(sympy.S(5)/2)
    F = -2*a**3/(3*x**(sympy.S(3)/2)) + 6*a**2*b*sqrt(x) + 6*a*b**2*x**(sympy.S(5)/2)/5 + 2*b**3*x**(sympy.S(9)/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_287():
    f = (a + b*x**2)**3/x**(sympy.S(7)/2)
    F = -2*a**3/(5*x**(sympy.S(5)/2)) - 6*a**2*b/sqrt(x) + 2*a*b**2*x**(sympy.S(3)/2) + 2*b**3*x**(sympy.S(7)/2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_288():
    f = x**(sympy.S(7)/2)/(a + b*x**2)
    F = -sqrt(2)*a**(sympy.S(5)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(9)/4)) + sqrt(2)*a**(sympy.S(5)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(9)/4)) - sqrt(2)*a**(sympy.S(5)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(9)/4)) + sqrt(2)*a**(sympy.S(5)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(9)/4)) - 2*a*sqrt(x)/b**2 + 2*x**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_289():
    f = x**(sympy.S(5)/2)/(a + b*x**2)
    F = -sqrt(2)*a**(sympy.S(3)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(7)/4)) + sqrt(2)*a**(sympy.S(3)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(7)/4)) + sqrt(2)*a**(sympy.S(3)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(7)/4)) - sqrt(2)*a**(sympy.S(3)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(7)/4)) + 2*x**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_290():
    f = x**(sympy.S(3)/2)/(a + b*x**2)
    F = sqrt(2)*a**(sympy.S(1)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(5)/4)) - sqrt(2)*a**(sympy.S(1)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(5)/4)) + sqrt(2)*a**(sympy.S(1)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(5)/4)) - sqrt(2)*a**(sympy.S(1)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(5)/4)) + 2*sqrt(x)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_291():
    f = sqrt(x)/(a + b*x**2)
    F = sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)) - sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)) - sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)) + sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_292():
    f = 1/(sqrt(x)*(a + b*x**2))
    F = -sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)) + sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)) - sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)) + sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_293():
    f = 1/(x**(sympy.S(3)/2)*(a + b*x**2))
    F = -2/(a*sqrt(x)) - sqrt(2)*b**(sympy.S(1)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(5)/4)) + sqrt(2)*b**(sympy.S(1)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(5)/4)) + sqrt(2)*b**(sympy.S(1)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)) - sqrt(2)*b**(sympy.S(1)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_294():
    f = 1/(x**(sympy.S(5)/2)*(a + b*x**2))
    F = -2/(3*a*x**(sympy.S(3)/2)) + sqrt(2)*b**(sympy.S(3)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(7)/4)) - sqrt(2)*b**(sympy.S(3)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(7)/4)) + sqrt(2)*b**(sympy.S(3)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4)) - sqrt(2)*b**(sympy.S(3)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_295():
    f = 1/(x**(sympy.S(7)/2)*(a + b*x**2))
    F = -2/(5*a*x**(sympy.S(5)/2)) + 2*b/(a**2*sqrt(x)) + sqrt(2)*b**(sympy.S(5)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(9)/4)) - sqrt(2)*b**(sympy.S(5)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(9)/4)) - sqrt(2)*b**(sympy.S(5)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(9)/4)) + sqrt(2)*b**(sympy.S(5)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_296():
    f = x**(sympy.S(7)/2)/(a + b*x**2)**2
    F = 5*sqrt(2)*a**(sympy.S(1)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*b**(sympy.S(9)/4)) - 5*sqrt(2)*a**(sympy.S(1)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*b**(sympy.S(9)/4)) + 5*sqrt(2)*a**(sympy.S(1)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*b**(sympy.S(9)/4)) - 5*sqrt(2)*a**(sympy.S(1)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*b**(sympy.S(9)/4)) - x**(sympy.S(5)/2)/(2*b*(a + b*x**2)) + 5*sqrt(x)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_297():
    f = x**(sympy.S(5)/2)/(a + b*x**2)**2
    F = -x**(sympy.S(3)/2)/(2*b*(a + b*x**2)) + 3*sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)) - 3*sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)) - 3*sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)) + 3*sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(1)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_298():
    f = x**(sympy.S(3)/2)/(a + b*x**2)**2
    F = -sqrt(x)/(2*b*(a + b*x**2)) - sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) + sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) - sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) + sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(3)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_299():
    f = sqrt(x)/(a + b*x**2)**2
    F = x**(sympy.S(3)/2)/(2*a*(a + b*x**2)) + sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(5)/4)*b**(sympy.S(3)/4)) - sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(5)/4)*b**(sympy.S(3)/4)) - sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*b**(sympy.S(3)/4)) + sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_300():
    f = 1/(sqrt(x)*(a + b*x**2)**2)
    F = sqrt(x)/(2*a*(a + b*x**2)) - 3*sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)) + 3*sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)) - 3*sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)) + 3*sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_301():
    f = 1/(x**(sympy.S(3)/2)*(a + b*x**2)**2)
    F = 1/(2*a*sqrt(x)*(a + b*x**2)) - 5/(2*a**2*sqrt(x)) - 5*sqrt(2)*b**(sympy.S(1)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(9)/4)) + 5*sqrt(2)*b**(sympy.S(1)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(9)/4)) + 5*sqrt(2)*b**(sympy.S(1)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(9)/4)) - 5*sqrt(2)*b**(sympy.S(1)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_302():
    f = 1/(x**(sympy.S(5)/2)*(a + b*x**2)**2)
    F = 1/(2*a*x**(sympy.S(3)/2)*(a + b*x**2)) - 7/(6*a**2*x**(sympy.S(3)/2)) + 7*sqrt(2)*b**(sympy.S(3)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(11)/4)) - 7*sqrt(2)*b**(sympy.S(3)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(11)/4)) + 7*sqrt(2)*b**(sympy.S(3)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(11)/4)) - 7*sqrt(2)*b**(sympy.S(3)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_303():
    f = 1/(x**(sympy.S(7)/2)*(a + b*x**2)**2)
    F = 1/(2*a*x**(sympy.S(5)/2)*(a + b*x**2)) - 9/(10*a**2*x**(sympy.S(5)/2)) + 9*b/(2*a**3*sqrt(x)) + 9*sqrt(2)*b**(sympy.S(5)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(13)/4)) - 9*sqrt(2)*b**(sympy.S(5)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(13)/4)) - 9*sqrt(2)*b**(sympy.S(5)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(13)/4)) + 9*sqrt(2)*b**(sympy.S(5)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_304():
    f = x**(sympy.S(7)/2)/(a + b*x**2)**3
    F = -x**(sympy.S(5)/2)/(4*b*(a + b*x**2)**2) - 5*sqrt(x)/(16*b**2*(a + b*x**2)) - 5*sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(3)/4)*b**(sympy.S(9)/4)) + 5*sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(3)/4)*b**(sympy.S(9)/4)) - 5*sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(3)/4)*b**(sympy.S(9)/4)) + 5*sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(3)/4)*b**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_305():
    f = x**(sympy.S(5)/2)/(a + b*x**2)**3
    F = -x**(sympy.S(3)/2)/(4*b*(a + b*x**2)**2) + 3*x**(sympy.S(3)/2)/(16*a*b*(a + b*x**2)) + 3*sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(5)/4)*b**(sympy.S(7)/4)) - 3*sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(5)/4)*b**(sympy.S(7)/4)) - 3*sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(5)/4)*b**(sympy.S(7)/4)) + 3*sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(5)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_306():
    f = x**(sympy.S(3)/2)/(a + b*x**2)**3
    F = -sqrt(x)/(4*b*(a + b*x**2)**2) + sqrt(x)/(16*a*b*(a + b*x**2)) - 3*sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) + 3*sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) - 3*sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) + 3*sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(7)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_307():
    f = sqrt(x)/(a + b*x**2)**3
    F = x**(sympy.S(3)/2)/(4*a*(a + b*x**2)**2) + 5*x**(sympy.S(3)/2)/(16*a**2*(a + b*x**2)) + 5*sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(9)/4)*b**(sympy.S(3)/4)) - 5*sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(9)/4)*b**(sympy.S(3)/4)) - 5*sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(9)/4)*b**(sympy.S(3)/4)) + 5*sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(9)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_308():
    f = 1/(sqrt(x)*(a + b*x**2)**3)
    F = sqrt(x)/(4*a*(a + b*x**2)**2) + 7*sqrt(x)/(16*a**2*(a + b*x**2)) - 21*sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(11)/4)*b**(sympy.S(1)/4)) + 21*sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(11)/4)*b**(sympy.S(1)/4)) - 21*sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*b**(sympy.S(1)/4)) + 21*sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*b**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_309():
    f = 1/(x**(sympy.S(3)/2)*(a + b*x**2)**3)
    F = 1/(4*a*sqrt(x)*(a + b*x**2)**2) + 9/(16*a**2*sqrt(x)*(a + b*x**2)) - 45/(16*a**3*sqrt(x)) - 45*sqrt(2)*b**(sympy.S(1)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(13)/4)) + 45*sqrt(2)*b**(sympy.S(1)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(13)/4)) + 45*sqrt(2)*b**(sympy.S(1)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(13)/4)) - 45*sqrt(2)*b**(sympy.S(1)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_310():
    f = 1/(x**(sympy.S(5)/2)*(a + b*x**2)**3)
    F = 1/(4*a*x**(sympy.S(3)/2)*(a + b*x**2)**2) + 11/(16*a**2*x**(sympy.S(3)/2)*(a + b*x**2)) - 77/(48*a**3*x**(sympy.S(3)/2)) + 77*sqrt(2)*b**(sympy.S(3)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(15)/4)) - 77*sqrt(2)*b**(sympy.S(3)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(15)/4)) + 77*sqrt(2)*b**(sympy.S(3)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(15)/4)) - 77*sqrt(2)*b**(sympy.S(3)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(15)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_311():
    f = 1/(x**(sympy.S(7)/2)*(a + b*x**2)**3)
    F = 1/(4*a*x**(sympy.S(5)/2)*(a + b*x**2)**2) + 13/(16*a**2*x**(sympy.S(5)/2)*(a + b*x**2)) - 117/(80*a**3*x**(sympy.S(5)/2)) + 117*b/(16*a**4*sqrt(x)) + 117*sqrt(2)*b**(sympy.S(5)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(17)/4)) - 117*sqrt(2)*b**(sympy.S(5)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(17)/4)) - 117*sqrt(2)*b**(sympy.S(5)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(17)/4)) + 117*sqrt(2)*b**(sympy.S(5)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(17)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_312():
    f = sqrt(x)/(a - b*x**2)
    F = -atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(a**(sympy.S(1)/4)*b**(sympy.S(3)/4)) + atanh(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(a**(sympy.S(1)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_313():
    f = x**(sympy.S(7)/2)/(x**2 + 1)
    F = 2*x**(sympy.S(5)/2)/5 - 2*sqrt(x) - sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/4 + sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/4 + sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/2 + sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_314():
    f = x**(sympy.S(5)/2)/(x**2 + 1)
    F = 2*x**(sympy.S(3)/2)/3 - sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/4 + sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/4 - sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/2 - sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_315():
    f = x**(sympy.S(3)/2)/(x**2 + 1)
    F = 2*sqrt(x) + sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/4 - sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/4 - sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/2 - sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_316():
    f = sqrt(x)/(x**2 + 1)
    F = sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/4 - sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/4 + sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/2 + sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_317():
    f = 1/(sqrt(x)*(x**2 + 1))
    F = -sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/4 + sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/4 + sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/2 + sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_318():
    f = 1/(x**(sympy.S(3)/2)*(x**2 + 1))
    F = -sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/4 + sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/4 - sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/2 - sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/2 - 2/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_319():
    f = 1/(x**(sympy.S(5)/2)*(x**2 + 1))
    F = sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/4 - sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/4 - sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/2 - sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/2 - 2/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_320():
    f = 1/(x**(sympy.S(7)/2)*(x**2 + 1))
    F = sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/4 - sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/4 + sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/2 + sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/2 + 2/sqrt(x) - 2/(5*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_321():
    f = x**(sympy.S(7)/2)/(x**2 + 1)**2
    F = -x**(sympy.S(5)/2)/(2*x**2 + 2) + 5*sqrt(x)/2 + 5*sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/16 - 5*sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/16 - 5*sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/8 - 5*sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_322():
    f = x**(sympy.S(5)/2)/(x**2 + 1)**2
    F = -x**(sympy.S(3)/2)/(2*x**2 + 2) + 3*sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/16 - 3*sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/16 + 3*sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/8 + 3*sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_323():
    f = x**(sympy.S(3)/2)/(x**2 + 1)**2
    F = -sqrt(x)/(2*x**2 + 2) - sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/16 + sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/16 + sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/8 + sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_324():
    f = sqrt(x)/(x**2 + 1)**2
    F = x**(sympy.S(3)/2)/(2*x**2 + 2) + sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/16 - sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/16 + sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/8 + sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_325():
    f = 1/(sqrt(x)*(x**2 + 1)**2)
    F = sqrt(x)/(2*x**2 + 2) - 3*sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/16 + 3*sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/16 + 3*sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/8 + 3*sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_326():
    f = 1/(x**(sympy.S(3)/2)*(x**2 + 1)**2)
    F = -5*sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/16 + 5*sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/16 - 5*sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/8 - 5*sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/8 - 5/(2*sqrt(x)) + 1/(2*sqrt(x)*(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_327():
    f = 1/(x**(sympy.S(5)/2)*(x**2 + 1)**2)
    F = 7*sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/16 - 7*sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/16 - 7*sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/8 - 7*sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/8 - 7/(6*x**(sympy.S(3)/2)) + 1/(2*x**(sympy.S(3)/2)*(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_328():
    f = 1/(x**(sympy.S(7)/2)*(x**2 + 1)**2)
    F = 9*sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/16 - 9*sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/16 + 9*sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/8 + 9*sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/8 + 9/(2*sqrt(x)) - 9/(10*x**(sympy.S(5)/2)) + 1/(2*x**(sympy.S(5)/2)*(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_329():
    f = x**(sympy.S(7)/2)/(x**2 + 1)**3
    F = -x**(sympy.S(5)/2)/(4*(x**2 + 1)**2) - 5*sqrt(x)/(16*x**2 + 16) - 5*sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/128 + 5*sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/128 + 5*sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/64 + 5*sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_330():
    f = x**(sympy.S(5)/2)/(x**2 + 1)**3
    F = 3*x**(sympy.S(3)/2)/(16*x**2 + 16) - x**(sympy.S(3)/2)/(4*(x**2 + 1)**2) + 3*sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/128 - 3*sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/128 + 3*sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/64 + 3*sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_331():
    f = x**(sympy.S(3)/2)/(x**2 + 1)**3
    F = sqrt(x)/(16*x**2 + 16) - sqrt(x)/(4*(x**2 + 1)**2) - 3*sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/128 + 3*sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/128 + 3*sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/64 + 3*sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_332():
    f = sqrt(x)/(x**2 + 1)**3
    F = 5*x**(sympy.S(3)/2)/(16*x**2 + 16) + x**(sympy.S(3)/2)/(4*(x**2 + 1)**2) + 5*sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/128 - 5*sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/128 + 5*sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/64 + 5*sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_333():
    f = 1/(sqrt(x)*(x**2 + 1)**3)
    F = 7*sqrt(x)/(16*x**2 + 16) + sqrt(x)/(4*(x**2 + 1)**2) - 21*sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/128 + 21*sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/128 + 21*sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/64 + 21*sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_334():
    f = 1/(x**(sympy.S(3)/2)*(x**2 + 1)**3)
    F = -45*sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/128 + 45*sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/128 - 45*sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/64 - 45*sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/64 - 45/(16*sqrt(x)) + 9/(16*sqrt(x)*(x**2 + 1)) + 1/(4*sqrt(x)*(x**2 + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_335():
    f = 1/(x**(sympy.S(5)/2)*(x**2 + 1)**3)
    F = 77*sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/128 - 77*sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/128 - 77*sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/64 - 77*sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/64 - 77/(48*x**(sympy.S(3)/2)) + 11/(16*x**(sympy.S(3)/2)*(x**2 + 1)) + 1/(4*x**(sympy.S(3)/2)*(x**2 + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_336():
    f = 1/(x**(sympy.S(7)/2)*(x**2 + 1)**3)
    F = 117*sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/128 - 117*sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/128 + 117*sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/64 + 117*sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/64 + 117/(16*sqrt(x)) - 117/(80*x**(sympy.S(5)/2)) + 13/(16*x**(sympy.S(5)/2)*(x**2 + 1)) + 1/(4*x**(sympy.S(5)/2)*(x**2 + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_337():
    f = sqrt(x)/(1 - x**2)
    F = -atan(sqrt(x)) + atanh(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_338():
    f = x**m*(a + b*x**2)**5
    F = a**5*x**(m + 1)/(m + 1) + 5*a**4*b*x**(m + 3)/(m + 3) + 10*a**3*b**2*x**(m + 5)/(m + 5) + 10*a**2*b**3*x**(m + 7)/(m + 7) + 5*a*b**4*x**(m + 9)/(m + 9) + b**5*x**(m + 11)/(m + 11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_339():
    f = x**m*(a + b*x**2)**4
    F = a**4*x**(m + 1)/(m + 1) + 4*a**3*b*x**(m + 3)/(m + 3) + 6*a**2*b**2*x**(m + 5)/(m + 5) + 4*a*b**3*x**(m + 7)/(m + 7) + b**4*x**(m + 9)/(m + 9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_340():
    f = x**m*(a + b*x**2)**3
    F = a**3*x**(m + 1)/(m + 1) + 3*a**2*b*x**(m + 3)/(m + 3) + 3*a*b**2*x**(m + 5)/(m + 5) + b**3*x**(m + 7)/(m + 7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_341():
    f = x**m*(a + b*x**2)**2
    F = a**2*x**(m + 1)/(m + 1) + 2*a*b*x**(m + 3)/(m + 3) + b**2*x**(m + 5)/(m + 5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_342():
    f = x**m*(a + b*x**2)
    F = a*x**(m + 1)/(m + 1) + b*x**(m + 3)/(m + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_343():
    f = x**m/(a + b*x**2)
    F = x**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_344():
    f = x**m/(a + b*x**2)**2
    F = x**(m + 1)*hyper((2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_345():
    f = x**m/(a + b*x**2)**3
    F = x**(m + 1)*hyper((3, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a**3*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_346():
    f = (c*x)**(m + 1)/(a + b*x**2)
    F = (c*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -b*x**2/a)/(a*c*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_347():
    f = (c*x)**m/(a + b*x**2)
    F = (c*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_348():
    f = (c*x)**(m - 1)/(a + b*x**2)
    F = (c*x)**m*hyper((1, m/2), (m/2 + 1,), -b*x**2/a)/(a*c*m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_349():
    f = (c*x)**(m - 2)/(a + b*x**2)
    F = -(c*x)**(m - 1)*hyper((1, m/2 + sympy.S(-1)/2), (m/2 + sympy.S.Half,), -b*x**2/a)/(a*c*(1 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_350():
    f = (c*x)**(m - 3)/(a + b*x**2)
    F = -(c*x)**(m - 2)*hyper((1, m/2 - 1), (m/2,), -b*x**2/a)/(a*c*(2 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_351():
    f = x**m/(a*x**2/b + 1)**2
    F = x**(m + 1)*hyper((2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -a*x**2/b)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_352():
    f = x**7*sqrt(a + b*x**2)
    F = -a**3*(a + b*x**2)**(sympy.S(3)/2)/(3*b**4) + 3*a**2*(a + b*x**2)**(sympy.S(5)/2)/(5*b**4) - 3*a*(a + b*x**2)**(sympy.S(7)/2)/(7*b**4) + (a + b*x**2)**(sympy.S(9)/2)/(9*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_353():
    f = x**5*sqrt(a + b*x**2)
    F = a**2*(a + b*x**2)**(sympy.S(3)/2)/(3*b**3) - 2*a*(a + b*x**2)**(sympy.S(5)/2)/(5*b**3) + (a + b*x**2)**(sympy.S(7)/2)/(7*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_354():
    f = x**3*sqrt(a + b*x**2)
    F = -a*(a + b*x**2)**(sympy.S(3)/2)/(3*b**2) + (a + b*x**2)**(sympy.S(5)/2)/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_355():
    f = x*sqrt(a + b*x**2)
    F = (a + b*x**2)**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_356():
    f = sqrt(a + b*x**2)/x
    F = -sqrt(a)*atanh(sqrt(a + b*x**2)/sqrt(a)) + sqrt(a + b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_357():
    f = sqrt(a + b*x**2)/x**3
    F = -sqrt(a + b*x**2)/(2*x**2) - b*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_358():
    f = sqrt(a + b*x**2)/x**5
    F = -sqrt(a + b*x**2)/(4*x**4) - b*sqrt(a + b*x**2)/(8*a*x**2) + b**2*atanh(sqrt(a + b*x**2)/sqrt(a))/(8*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_359():
    f = sqrt(a + b*x**2)/x**7
    F = -sqrt(a + b*x**2)/(6*x**6) - b*sqrt(a + b*x**2)/(24*a*x**4) + b**2*sqrt(a + b*x**2)/(16*a**2*x**2) - b**3*atanh(sqrt(a + b*x**2)/sqrt(a))/(16*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_360():
    f = x**4*sqrt(a + b*x**2)
    F = a**3*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*b**(sympy.S(5)/2)) - a**2*x*sqrt(a + b*x**2)/(16*b**2) + a*x**3*sqrt(a + b*x**2)/(24*b) + x**5*sqrt(a + b*x**2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_361():
    f = x**2*sqrt(a + b*x**2)
    F = -a**2*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(3)/2)) + a*x*sqrt(a + b*x**2)/(8*b) + x**3*sqrt(a + b*x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_362():
    f = sqrt(a + b*x**2)
    F = a*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*sqrt(b)) + x*sqrt(a + b*x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_363():
    f = sqrt(a + b*x**2)/x**2
    F = sqrt(b)*atanh(sqrt(b)*x/sqrt(a + b*x**2)) - sqrt(a + b*x**2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_364():
    f = sqrt(a + b*x**2)/x**4
    F = -(a + b*x**2)**(sympy.S(3)/2)/(3*a*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_365():
    f = sqrt(a + b*x**2)/x**6
    F = -(a + b*x**2)**(sympy.S(3)/2)/(5*a*x**5) + 2*b*(a + b*x**2)**(sympy.S(3)/2)/(15*a**2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_366():
    f = sqrt(a + b*x**2)/x**8
    F = -(a + b*x**2)**(sympy.S(3)/2)/(7*a*x**7) + 4*b*(a + b*x**2)**(sympy.S(3)/2)/(35*a**2*x**5) - 8*b**2*(a + b*x**2)**(sympy.S(3)/2)/(105*a**3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_367():
    f = sqrt(a + b*x**2)/x**10
    F = -(a + b*x**2)**(sympy.S(3)/2)/(9*a*x**9) + 2*b*(a + b*x**2)**(sympy.S(3)/2)/(21*a**2*x**7) - 8*b**2*(a + b*x**2)**(sympy.S(3)/2)/(105*a**3*x**5) + 16*b**3*(a + b*x**2)**(sympy.S(3)/2)/(315*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_368():
    f = x**7*(a + b*x**2)**(sympy.S(3)/2)
    F = -a**3*(a + b*x**2)**(sympy.S(5)/2)/(5*b**4) + 3*a**2*(a + b*x**2)**(sympy.S(7)/2)/(7*b**4) - a*(a + b*x**2)**(sympy.S(9)/2)/(3*b**4) + (a + b*x**2)**(sympy.S(11)/2)/(11*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_369():
    f = x**5*(a + b*x**2)**(sympy.S(3)/2)
    F = a**2*(a + b*x**2)**(sympy.S(5)/2)/(5*b**3) - 2*a*(a + b*x**2)**(sympy.S(7)/2)/(7*b**3) + (a + b*x**2)**(sympy.S(9)/2)/(9*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_370():
    f = x**3*(a + b*x**2)**(sympy.S(3)/2)
    F = -a*(a + b*x**2)**(sympy.S(5)/2)/(5*b**2) + (a + b*x**2)**(sympy.S(7)/2)/(7*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_371():
    f = x*(a + b*x**2)**(sympy.S(3)/2)
    F = (a + b*x**2)**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_372():
    f = (a + b*x**2)**(sympy.S(3)/2)/x
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a + b*x**2)/sqrt(a)) + a*sqrt(a + b*x**2) + (a + b*x**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_373():
    f = (a + b*x**2)**(sympy.S(3)/2)/x**3
    F = -3*sqrt(a)*b*atanh(sqrt(a + b*x**2)/sqrt(a))/2 + 3*b*sqrt(a + b*x**2)/2 - (a + b*x**2)**(sympy.S(3)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_374():
    f = (a + b*x**2)**(sympy.S(3)/2)/x**5
    F = -3*b*sqrt(a + b*x**2)/(8*x**2) - (a + b*x**2)**(sympy.S(3)/2)/(4*x**4) - 3*b**2*atanh(sqrt(a + b*x**2)/sqrt(a))/(8*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_375():
    f = (a + b*x**2)**(sympy.S(3)/2)/x**7
    F = -b*sqrt(a + b*x**2)/(8*x**4) - (a + b*x**2)**(sympy.S(3)/2)/(6*x**6) - b**2*sqrt(a + b*x**2)/(16*a*x**2) + b**3*atanh(sqrt(a + b*x**2)/sqrt(a))/(16*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_376():
    f = (a + b*x**2)**(sympy.S(3)/2)/x**9
    F = -b*sqrt(a + b*x**2)/(16*x**6) - (a + b*x**2)**(sympy.S(3)/2)/(8*x**8) - b**2*sqrt(a + b*x**2)/(64*a*x**4) + 3*b**3*sqrt(a + b*x**2)/(128*a**2*x**2) - 3*b**4*atanh(sqrt(a + b*x**2)/sqrt(a))/(128*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_377():
    f = x**4*(a + b*x**2)**(sympy.S(3)/2)
    F = 3*a**4*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(128*b**(sympy.S(5)/2)) - 3*a**3*x*sqrt(a + b*x**2)/(128*b**2) + a**2*x**3*sqrt(a + b*x**2)/(64*b) + a*x**5*sqrt(a + b*x**2)/16 + x**5*(a + b*x**2)**(sympy.S(3)/2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_378():
    f = x**2*(a + b*x**2)**(sympy.S(3)/2)
    F = -a**3*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*b**(sympy.S(3)/2)) + a**2*x*sqrt(a + b*x**2)/(16*b) + a*x**3*sqrt(a + b*x**2)/8 + x**3*(a + b*x**2)**(sympy.S(3)/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_379():
    f = (a + b*x**2)**(sympy.S(3)/2)
    F = 3*a**2*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*sqrt(b)) + 3*a*x*sqrt(a + b*x**2)/8 + x*(a + b*x**2)**(sympy.S(3)/2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_380():
    f = (a + b*x**2)**(sympy.S(3)/2)/x**2
    F = 3*a*sqrt(b)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/2 + 3*b*x*sqrt(a + b*x**2)/2 - (a + b*x**2)**(sympy.S(3)/2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_381():
    f = (a + b*x**2)**(sympy.S(3)/2)/x**4
    F = b**(sympy.S(3)/2)*atanh(sqrt(b)*x/sqrt(a + b*x**2)) - b*sqrt(a + b*x**2)/x - (a + b*x**2)**(sympy.S(3)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_382():
    f = (a + b*x**2)**(sympy.S(3)/2)/x**6
    F = -(a + b*x**2)**(sympy.S(5)/2)/(5*a*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_383():
    f = (a + b*x**2)**(sympy.S(3)/2)/x**8
    F = -(a + b*x**2)**(sympy.S(5)/2)/(7*a*x**7) + 2*b*(a + b*x**2)**(sympy.S(5)/2)/(35*a**2*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_384():
    f = (a + b*x**2)**(sympy.S(3)/2)/x**10
    F = -(a + b*x**2)**(sympy.S(5)/2)/(9*a*x**9) + 4*b*(a + b*x**2)**(sympy.S(5)/2)/(63*a**2*x**7) - 8*b**2*(a + b*x**2)**(sympy.S(5)/2)/(315*a**3*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_385():
    f = (a + b*x**2)**(sympy.S(3)/2)/x**12
    F = -(a + b*x**2)**(sympy.S(5)/2)/(11*a*x**11) + 2*b*(a + b*x**2)**(sympy.S(5)/2)/(33*a**2*x**9) - 8*b**2*(a + b*x**2)**(sympy.S(5)/2)/(231*a**3*x**7) + 16*b**3*(a + b*x**2)**(sympy.S(5)/2)/(1155*a**4*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_386():
    f = x**7*(a + b*x**2)**(sympy.S(5)/2)
    F = -a**3*(a + b*x**2)**(sympy.S(7)/2)/(7*b**4) + a**2*(a + b*x**2)**(sympy.S(9)/2)/(3*b**4) - 3*a*(a + b*x**2)**(sympy.S(11)/2)/(11*b**4) + (a + b*x**2)**(sympy.S(13)/2)/(13*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_387():
    f = x**5*(a + b*x**2)**(sympy.S(5)/2)
    F = a**2*(a + b*x**2)**(sympy.S(7)/2)/(7*b**3) - 2*a*(a + b*x**2)**(sympy.S(9)/2)/(9*b**3) + (a + b*x**2)**(sympy.S(11)/2)/(11*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_388():
    f = x**3*(a + b*x**2)**(sympy.S(5)/2)
    F = -a*(a + b*x**2)**(sympy.S(7)/2)/(7*b**2) + (a + b*x**2)**(sympy.S(9)/2)/(9*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_389():
    f = x*(a + b*x**2)**(sympy.S(5)/2)
    F = (a + b*x**2)**(sympy.S(7)/2)/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_390():
    f = (a + b*x**2)**(sympy.S(5)/2)/x
    F = -a**(sympy.S(5)/2)*atanh(sqrt(a + b*x**2)/sqrt(a)) + a**2*sqrt(a + b*x**2) + a*(a + b*x**2)**(sympy.S(3)/2)/3 + (a + b*x**2)**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_391():
    f = (a + b*x**2)**(sympy.S(5)/2)/x**3
    F = -5*a**(sympy.S(3)/2)*b*atanh(sqrt(a + b*x**2)/sqrt(a))/2 + 5*a*b*sqrt(a + b*x**2)/2 + 5*b*(a + b*x**2)**(sympy.S(3)/2)/6 - (a + b*x**2)**(sympy.S(5)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_392():
    f = (a + b*x**2)**(sympy.S(5)/2)/x**5
    F = -15*sqrt(a)*b**2*atanh(sqrt(a + b*x**2)/sqrt(a))/8 + 15*b**2*sqrt(a + b*x**2)/8 - 5*b*(a + b*x**2)**(sympy.S(3)/2)/(8*x**2) - (a + b*x**2)**(sympy.S(5)/2)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_393():
    f = (a + b*x**2)**(sympy.S(5)/2)/x**7
    F = -5*b**2*sqrt(a + b*x**2)/(16*x**2) - 5*b*(a + b*x**2)**(sympy.S(3)/2)/(24*x**4) - (a + b*x**2)**(sympy.S(5)/2)/(6*x**6) - 5*b**3*atanh(sqrt(a + b*x**2)/sqrt(a))/(16*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_394():
    f = (a + b*x**2)**(sympy.S(5)/2)/x**9
    F = -5*b**2*sqrt(a + b*x**2)/(64*x**4) - 5*b*(a + b*x**2)**(sympy.S(3)/2)/(48*x**6) - (a + b*x**2)**(sympy.S(5)/2)/(8*x**8) - 5*b**3*sqrt(a + b*x**2)/(128*a*x**2) + 5*b**4*atanh(sqrt(a + b*x**2)/sqrt(a))/(128*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_395():
    f = (a + b*x**2)**(sympy.S(5)/2)/x**11
    F = -b**2*sqrt(a + b*x**2)/(32*x**6) - b*(a + b*x**2)**(sympy.S(3)/2)/(16*x**8) - (a + b*x**2)**(sympy.S(5)/2)/(10*x**10) - b**3*sqrt(a + b*x**2)/(128*a*x**4) + 3*b**4*sqrt(a + b*x**2)/(256*a**2*x**2) - 3*b**5*atanh(sqrt(a + b*x**2)/sqrt(a))/(256*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_396():
    f = x**4*(a + b*x**2)**(sympy.S(5)/2)
    F = 3*a**5*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(256*b**(sympy.S(5)/2)) - 3*a**4*x*sqrt(a + b*x**2)/(256*b**2) + a**3*x**3*sqrt(a + b*x**2)/(128*b) + a**2*x**5*sqrt(a + b*x**2)/32 + a*x**5*(a + b*x**2)**(sympy.S(3)/2)/16 + x**5*(a + b*x**2)**(sympy.S(5)/2)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_397():
    f = x**2*(a + b*x**2)**(sympy.S(5)/2)
    F = -5*a**4*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(128*b**(sympy.S(3)/2)) + 5*a**3*x*sqrt(a + b*x**2)/(128*b) + 5*a**2*x**3*sqrt(a + b*x**2)/64 + 5*a*x**3*(a + b*x**2)**(sympy.S(3)/2)/48 + x**3*(a + b*x**2)**(sympy.S(5)/2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_398():
    f = (a + b*x**2)**(sympy.S(5)/2)
    F = 5*a**3*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*sqrt(b)) + 5*a**2*x*sqrt(a + b*x**2)/16 + 5*a*x*(a + b*x**2)**(sympy.S(3)/2)/24 + x*(a + b*x**2)**(sympy.S(5)/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_399():
    f = (a + b*x**2)**(sympy.S(5)/2)/x**2
    F = 15*a**2*sqrt(b)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/8 + 15*a*b*x*sqrt(a + b*x**2)/8 + 5*b*x*(a + b*x**2)**(sympy.S(3)/2)/4 - (a + b*x**2)**(sympy.S(5)/2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_400():
    f = (a + b*x**2)**(sympy.S(5)/2)/x**4
    F = 5*a*b**(sympy.S(3)/2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/2 + 5*b**2*x*sqrt(a + b*x**2)/2 - 5*b*(a + b*x**2)**(sympy.S(3)/2)/(3*x) - (a + b*x**2)**(sympy.S(5)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_401():
    f = (a + b*x**2)**(sympy.S(5)/2)/x**6
    F = b**(sympy.S(5)/2)*atanh(sqrt(b)*x/sqrt(a + b*x**2)) - b**2*sqrt(a + b*x**2)/x - b*(a + b*x**2)**(sympy.S(3)/2)/(3*x**3) - (a + b*x**2)**(sympy.S(5)/2)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_402():
    f = (a + b*x**2)**(sympy.S(5)/2)/x**8
    F = -(a + b*x**2)**(sympy.S(7)/2)/(7*a*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_403():
    f = (a + b*x**2)**(sympy.S(5)/2)/x**10
    F = -(a + b*x**2)**(sympy.S(7)/2)/(9*a*x**9) + 2*b*(a + b*x**2)**(sympy.S(7)/2)/(63*a**2*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_404():
    f = (a + b*x**2)**(sympy.S(5)/2)/x**12
    F = -(a + b*x**2)**(sympy.S(7)/2)/(11*a*x**11) + 4*b*(a + b*x**2)**(sympy.S(7)/2)/(99*a**2*x**9) - 8*b**2*(a + b*x**2)**(sympy.S(7)/2)/(693*a**3*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_405():
    f = (a + b*x**2)**(sympy.S(5)/2)/x**14
    F = -(a + b*x**2)**(sympy.S(7)/2)/(13*a*x**13) + 6*b*(a + b*x**2)**(sympy.S(7)/2)/(143*a**2*x**11) - 8*b**2*(a + b*x**2)**(sympy.S(7)/2)/(429*a**3*x**9) + 16*b**3*(a + b*x**2)**(sympy.S(7)/2)/(3003*a**4*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_406():
    f = (a + b*x**2)**(sympy.S(5)/2)/x**16
    F = -(a + b*x**2)**(sympy.S(7)/2)/(15*a*x**15) + 8*b*(a + b*x**2)**(sympy.S(7)/2)/(195*a**2*x**13) - 16*b**2*(a + b*x**2)**(sympy.S(7)/2)/(715*a**3*x**11) + 64*b**3*(a + b*x**2)**(sympy.S(7)/2)/(6435*a**4*x**9) - 128*b**4*(a + b*x**2)**(sympy.S(7)/2)/(45045*a**5*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_407():
    f = (a + b*x**2)**(sympy.S(5)/2)/x**18
    F = -(a + b*x**2)**(sympy.S(7)/2)/(17*a*x**17) + 2*b*(a + b*x**2)**(sympy.S(7)/2)/(51*a**2*x**15) - 16*b**2*(a + b*x**2)**(sympy.S(7)/2)/(663*a**3*x**13) + 32*b**3*(a + b*x**2)**(sympy.S(7)/2)/(2431*a**4*x**11) - 128*b**4*(a + b*x**2)**(sympy.S(7)/2)/(21879*a**5*x**9) + 256*b**5*(a + b*x**2)**(sympy.S(7)/2)/(153153*a**6*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_408():
    f = x**15*(a + b*x**2)**(sympy.S(9)/2)
    F = -a**7*(a + b*x**2)**(sympy.S(11)/2)/(11*b**8) + 7*a**6*(a + b*x**2)**(sympy.S(13)/2)/(13*b**8) - 7*a**5*(a + b*x**2)**(sympy.S(15)/2)/(5*b**8) + 35*a**4*(a + b*x**2)**(sympy.S(17)/2)/(17*b**8) - 35*a**3*(a + b*x**2)**(sympy.S(19)/2)/(19*b**8) + a**2*(a + b*x**2)**(sympy.S(21)/2)/b**8 - 7*a*(a + b*x**2)**(sympy.S(23)/2)/(23*b**8) + (a + b*x**2)**(sympy.S(25)/2)/(25*b**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_409():
    f = x**13*(a + b*x**2)**(sympy.S(9)/2)
    F = a**6*(a + b*x**2)**(sympy.S(11)/2)/(11*b**7) - 6*a**5*(a + b*x**2)**(sympy.S(13)/2)/(13*b**7) + a**4*(a + b*x**2)**(sympy.S(15)/2)/b**7 - 20*a**3*(a + b*x**2)**(sympy.S(17)/2)/(17*b**7) + 15*a**2*(a + b*x**2)**(sympy.S(19)/2)/(19*b**7) - 2*a*(a + b*x**2)**(sympy.S(21)/2)/(7*b**7) + (a + b*x**2)**(sympy.S(23)/2)/(23*b**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_410():
    f = x**11*(a + b*x**2)**(sympy.S(9)/2)
    F = -a**5*(a + b*x**2)**(sympy.S(11)/2)/(11*b**6) + 5*a**4*(a + b*x**2)**(sympy.S(13)/2)/(13*b**6) - 2*a**3*(a + b*x**2)**(sympy.S(15)/2)/(3*b**6) + 10*a**2*(a + b*x**2)**(sympy.S(17)/2)/(17*b**6) - 5*a*(a + b*x**2)**(sympy.S(19)/2)/(19*b**6) + (a + b*x**2)**(sympy.S(21)/2)/(21*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_411():
    f = x**9*(a + b*x**2)**(sympy.S(9)/2)
    F = a**4*(a + b*x**2)**(sympy.S(11)/2)/(11*b**5) - 4*a**3*(a + b*x**2)**(sympy.S(13)/2)/(13*b**5) + 2*a**2*(a + b*x**2)**(sympy.S(15)/2)/(5*b**5) - 4*a*(a + b*x**2)**(sympy.S(17)/2)/(17*b**5) + (a + b*x**2)**(sympy.S(19)/2)/(19*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_412():
    f = x**7*(a + b*x**2)**(sympy.S(9)/2)
    F = -a**3*(a + b*x**2)**(sympy.S(11)/2)/(11*b**4) + 3*a**2*(a + b*x**2)**(sympy.S(13)/2)/(13*b**4) - a*(a + b*x**2)**(sympy.S(15)/2)/(5*b**4) + (a + b*x**2)**(sympy.S(17)/2)/(17*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_413():
    f = x**5*(a + b*x**2)**(sympy.S(9)/2)
    F = a**2*(a + b*x**2)**(sympy.S(11)/2)/(11*b**3) - 2*a*(a + b*x**2)**(sympy.S(13)/2)/(13*b**3) + (a + b*x**2)**(sympy.S(15)/2)/(15*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_414():
    f = x**3*(a + b*x**2)**(sympy.S(9)/2)
    F = -a*(a + b*x**2)**(sympy.S(11)/2)/(11*b**2) + (a + b*x**2)**(sympy.S(13)/2)/(13*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_415():
    f = x*(a + b*x**2)**(sympy.S(9)/2)
    F = (a + b*x**2)**(sympy.S(11)/2)/(11*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_416():
    f = (a + b*x**2)**(sympy.S(9)/2)/x
    F = -a**(sympy.S(9)/2)*atanh(sqrt(a + b*x**2)/sqrt(a)) + a**4*sqrt(a + b*x**2) + a**3*(a + b*x**2)**(sympy.S(3)/2)/3 + a**2*(a + b*x**2)**(sympy.S(5)/2)/5 + a*(a + b*x**2)**(sympy.S(7)/2)/7 + (a + b*x**2)**(sympy.S(9)/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_417():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**3
    F = -9*a**(sympy.S(7)/2)*b*atanh(sqrt(a + b*x**2)/sqrt(a))/2 + 9*a**3*b*sqrt(a + b*x**2)/2 + 3*a**2*b*(a + b*x**2)**(sympy.S(3)/2)/2 + 9*a*b*(a + b*x**2)**(sympy.S(5)/2)/10 + 9*b*(a + b*x**2)**(sympy.S(7)/2)/14 - (a + b*x**2)**(sympy.S(9)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_418():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**5
    F = -63*a**(sympy.S(5)/2)*b**2*atanh(sqrt(a + b*x**2)/sqrt(a))/8 + 63*a**2*b**2*sqrt(a + b*x**2)/8 + 21*a*b**2*(a + b*x**2)**(sympy.S(3)/2)/8 + 63*b**2*(a + b*x**2)**(sympy.S(5)/2)/40 - 9*b*(a + b*x**2)**(sympy.S(7)/2)/(8*x**2) - (a + b*x**2)**(sympy.S(9)/2)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_419():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**7
    F = -105*a**(sympy.S(3)/2)*b**3*atanh(sqrt(a + b*x**2)/sqrt(a))/16 + 105*a*b**3*sqrt(a + b*x**2)/16 + 35*b**3*(a + b*x**2)**(sympy.S(3)/2)/16 - 21*b**2*(a + b*x**2)**(sympy.S(5)/2)/(16*x**2) - 3*b*(a + b*x**2)**(sympy.S(7)/2)/(8*x**4) - (a + b*x**2)**(sympy.S(9)/2)/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_420():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**9
    F = -315*sqrt(a)*b**4*atanh(sqrt(a + b*x**2)/sqrt(a))/128 + 315*b**4*sqrt(a + b*x**2)/128 - 105*b**3*(a + b*x**2)**(sympy.S(3)/2)/(128*x**2) - 21*b**2*(a + b*x**2)**(sympy.S(5)/2)/(64*x**4) - 3*b*(a + b*x**2)**(sympy.S(7)/2)/(16*x**6) - (a + b*x**2)**(sympy.S(9)/2)/(8*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_421():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**11
    F = -63*b**4*sqrt(a + b*x**2)/(256*x**2) - 21*b**3*(a + b*x**2)**(sympy.S(3)/2)/(128*x**4) - 21*b**2*(a + b*x**2)**(sympy.S(5)/2)/(160*x**6) - 9*b*(a + b*x**2)**(sympy.S(7)/2)/(80*x**8) - (a + b*x**2)**(sympy.S(9)/2)/(10*x**10) - 63*b**5*atanh(sqrt(a + b*x**2)/sqrt(a))/(256*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_422():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**13
    F = -21*b**4*sqrt(a + b*x**2)/(512*x**4) - 7*b**3*(a + b*x**2)**(sympy.S(3)/2)/(128*x**6) - 21*b**2*(a + b*x**2)**(sympy.S(5)/2)/(320*x**8) - 3*b*(a + b*x**2)**(sympy.S(7)/2)/(40*x**10) - (a + b*x**2)**(sympy.S(9)/2)/(12*x**12) - 21*b**5*sqrt(a + b*x**2)/(1024*a*x**2) + 21*b**6*atanh(sqrt(a + b*x**2)/sqrt(a))/(1024*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_423():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**15
    F = -3*b**4*sqrt(a + b*x**2)/(256*x**6) - 3*b**3*(a + b*x**2)**(sympy.S(3)/2)/(128*x**8) - 3*b**2*(a + b*x**2)**(sympy.S(5)/2)/(80*x**10) - 3*b*(a + b*x**2)**(sympy.S(7)/2)/(56*x**12) - (a + b*x**2)**(sympy.S(9)/2)/(14*x**14) - 3*b**5*sqrt(a + b*x**2)/(1024*a*x**4) + 9*b**6*sqrt(a + b*x**2)/(2048*a**2*x**2) - 9*b**7*atanh(sqrt(a + b*x**2)/sqrt(a))/(2048*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_424():
    f = x**6*(a + b*x**2)**(sympy.S(9)/2)
    F = -45*a**8*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(32768*b**(sympy.S(7)/2)) + 45*a**7*x*sqrt(a + b*x**2)/(32768*b**3) - 15*a**6*x**3*sqrt(a + b*x**2)/(16384*b**2) + 3*a**5*x**5*sqrt(a + b*x**2)/(4096*b) + 9*a**4*x**7*sqrt(a + b*x**2)/2048 + 3*a**3*x**7*(a + b*x**2)**(sympy.S(3)/2)/256 + 3*a**2*x**7*(a + b*x**2)**(sympy.S(5)/2)/128 + 9*a*x**7*(a + b*x**2)**(sympy.S(7)/2)/224 + x**7*(a + b*x**2)**(sympy.S(9)/2)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_425():
    f = x**4*(a + b*x**2)**(sympy.S(9)/2)
    F = 9*a**7*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2048*b**(sympy.S(5)/2)) - 9*a**6*x*sqrt(a + b*x**2)/(2048*b**2) + 3*a**5*x**3*sqrt(a + b*x**2)/(1024*b) + 3*a**4*x**5*sqrt(a + b*x**2)/256 + 3*a**3*x**5*(a + b*x**2)**(sympy.S(3)/2)/128 + 3*a**2*x**5*(a + b*x**2)**(sympy.S(5)/2)/80 + 3*a*x**5*(a + b*x**2)**(sympy.S(7)/2)/56 + x**5*(a + b*x**2)**(sympy.S(9)/2)/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_426():
    f = x**2*(a + b*x**2)**(sympy.S(9)/2)
    F = -21*a**6*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(1024*b**(sympy.S(3)/2)) + 21*a**5*x*sqrt(a + b*x**2)/(1024*b) + 21*a**4*x**3*sqrt(a + b*x**2)/512 + 7*a**3*x**3*(a + b*x**2)**(sympy.S(3)/2)/128 + 21*a**2*x**3*(a + b*x**2)**(sympy.S(5)/2)/320 + 3*a*x**3*(a + b*x**2)**(sympy.S(7)/2)/40 + x**3*(a + b*x**2)**(sympy.S(9)/2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_427():
    f = (a + b*x**2)**(sympy.S(9)/2)
    F = 63*a**5*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(256*sqrt(b)) + 63*a**4*x*sqrt(a + b*x**2)/256 + 21*a**3*x*(a + b*x**2)**(sympy.S(3)/2)/128 + 21*a**2*x*(a + b*x**2)**(sympy.S(5)/2)/160 + 9*a*x*(a + b*x**2)**(sympy.S(7)/2)/80 + x*(a + b*x**2)**(sympy.S(9)/2)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_428():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**2
    F = 315*a**4*sqrt(b)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/128 + 315*a**3*b*x*sqrt(a + b*x**2)/128 + 105*a**2*b*x*(a + b*x**2)**(sympy.S(3)/2)/64 + 21*a*b*x*(a + b*x**2)**(sympy.S(5)/2)/16 + 9*b*x*(a + b*x**2)**(sympy.S(7)/2)/8 - (a + b*x**2)**(sympy.S(9)/2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_429():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**4
    F = 105*a**3*b**(sympy.S(3)/2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/16 + 105*a**2*b**2*x*sqrt(a + b*x**2)/16 + 35*a*b**2*x*(a + b*x**2)**(sympy.S(3)/2)/8 + 7*b**2*x*(a + b*x**2)**(sympy.S(5)/2)/2 - 3*b*(a + b*x**2)**(sympy.S(7)/2)/x - (a + b*x**2)**(sympy.S(9)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_430():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**6
    F = 63*a**2*b**(sympy.S(5)/2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/8 + 63*a*b**3*x*sqrt(a + b*x**2)/8 + 21*b**3*x*(a + b*x**2)**(sympy.S(3)/2)/4 - 21*b**2*(a + b*x**2)**(sympy.S(5)/2)/(5*x) - 3*b*(a + b*x**2)**(sympy.S(7)/2)/(5*x**3) - (a + b*x**2)**(sympy.S(9)/2)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_431():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**8
    F = 9*a*b**(sympy.S(7)/2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/2 + 9*b**4*x*sqrt(a + b*x**2)/2 - 3*b**3*(a + b*x**2)**(sympy.S(3)/2)/x - 3*b**2*(a + b*x**2)**(sympy.S(5)/2)/(5*x**3) - 9*b*(a + b*x**2)**(sympy.S(7)/2)/(35*x**5) - (a + b*x**2)**(sympy.S(9)/2)/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_432():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**10
    F = b**(sympy.S(9)/2)*atanh(sqrt(b)*x/sqrt(a + b*x**2)) - b**4*sqrt(a + b*x**2)/x - b**3*(a + b*x**2)**(sympy.S(3)/2)/(3*x**3) - b**2*(a + b*x**2)**(sympy.S(5)/2)/(5*x**5) - b*(a + b*x**2)**(sympy.S(7)/2)/(7*x**7) - (a + b*x**2)**(sympy.S(9)/2)/(9*x**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_433():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**12
    F = -(a + b*x**2)**(sympy.S(11)/2)/(11*a*x**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_434():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**14
    F = -(a + b*x**2)**(sympy.S(11)/2)/(13*a*x**13) + 2*b*(a + b*x**2)**(sympy.S(11)/2)/(143*a**2*x**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_435():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**16
    F = -(a + b*x**2)**(sympy.S(11)/2)/(15*a*x**15) + 4*b*(a + b*x**2)**(sympy.S(11)/2)/(195*a**2*x**13) - 8*b**2*(a + b*x**2)**(sympy.S(11)/2)/(2145*a**3*x**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_436():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**18
    F = -(a + b*x**2)**(sympy.S(11)/2)/(17*a*x**17) + 2*b*(a + b*x**2)**(sympy.S(11)/2)/(85*a**2*x**15) - 8*b**2*(a + b*x**2)**(sympy.S(11)/2)/(1105*a**3*x**13) + 16*b**3*(a + b*x**2)**(sympy.S(11)/2)/(12155*a**4*x**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_437():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**20
    F = -(a + b*x**2)**(sympy.S(11)/2)/(19*a*x**19) + 8*b*(a + b*x**2)**(sympy.S(11)/2)/(323*a**2*x**17) - 16*b**2*(a + b*x**2)**(sympy.S(11)/2)/(1615*a**3*x**15) + 64*b**3*(a + b*x**2)**(sympy.S(11)/2)/(20995*a**4*x**13) - 128*b**4*(a + b*x**2)**(sympy.S(11)/2)/(230945*a**5*x**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_438():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**22
    F = -(a + b*x**2)**(sympy.S(11)/2)/(21*a*x**21) + 10*b*(a + b*x**2)**(sympy.S(11)/2)/(399*a**2*x**19) - 80*b**2*(a + b*x**2)**(sympy.S(11)/2)/(6783*a**3*x**17) + 32*b**3*(a + b*x**2)**(sympy.S(11)/2)/(6783*a**4*x**15) - 128*b**4*(a + b*x**2)**(sympy.S(11)/2)/(88179*a**5*x**13) + 256*b**5*(a + b*x**2)**(sympy.S(11)/2)/(969969*a**6*x**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_439():
    f = (a + b*x**2)**(sympy.S(9)/2)/x**24
    F = -(a + b*x**2)**(sympy.S(11)/2)/(23*a*x**23) + 4*b*(a + b*x**2)**(sympy.S(11)/2)/(161*a**2*x**21) - 40*b**2*(a + b*x**2)**(sympy.S(11)/2)/(3059*a**3*x**19) + 320*b**3*(a + b*x**2)**(sympy.S(11)/2)/(52003*a**4*x**17) - 128*b**4*(a + b*x**2)**(sympy.S(11)/2)/(52003*a**5*x**15) + 512*b**5*(a + b*x**2)**(sympy.S(11)/2)/(676039*a**6*x**13) - 1024*b**6*(a + b*x**2)**(sympy.S(11)/2)/(7436429*a**7*x**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_440():
    f = x**5*sqrt(4*x**2 + 9)
    F = (4*x**2 + 9)**(sympy.S(7)/2)/448 - 9*(4*x**2 + 9)**(sympy.S(5)/2)/160 + 27*(4*x**2 + 9)**(sympy.S(3)/2)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_441():
    f = x**4*sqrt(4*x**2 + 9)
    F = x**5*sqrt(4*x**2 + 9)/6 + 3*x**3*sqrt(4*x**2 + 9)/32 - 81*x*sqrt(4*x**2 + 9)/256 + 729*asinh(2*x/3)/512
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_442():
    f = x**3*sqrt(4*x**2 + 9)
    F = (4*x**2 + 9)**(sympy.S(5)/2)/80 - 3*(4*x**2 + 9)**(sympy.S(3)/2)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_443():
    f = x**2*sqrt(4*x**2 + 9)
    F = x**3*sqrt(4*x**2 + 9)/4 + 9*x*sqrt(4*x**2 + 9)/32 - 81*asinh(2*x/3)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_444():
    f = x*sqrt(4*x**2 + 9)
    F = (4*x**2 + 9)**(sympy.S(3)/2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_445():
    f = sqrt(4*x**2 + 9)
    F = x*sqrt(4*x**2 + 9)/2 + 9*asinh(2*x/3)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_446():
    f = sqrt(4*x**2 + 9)/x
    F = sqrt(4*x**2 + 9) - 3*atanh(sqrt(4*x**2 + 9)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_447():
    f = sqrt(4*x**2 + 9)/x**2
    F = 2*asinh(2*x/3) - sqrt(4*x**2 + 9)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_448():
    f = sqrt(4*x**2 + 9)/x**3
    F = -2*atanh(sqrt(4*x**2 + 9)/3)/3 - sqrt(4*x**2 + 9)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_449():
    f = sqrt(4*x**2 + 9)/x**4
    F = -(4*x**2 + 9)**(sympy.S(3)/2)/(27*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_450():
    f = sqrt(4*x**2 + 9)/x**5
    F = 2*atanh(sqrt(4*x**2 + 9)/3)/27 - sqrt(4*x**2 + 9)/(18*x**2) - sqrt(4*x**2 + 9)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_451():
    f = x**5*sqrt(9 - 4*x**2)
    F = -(9 - 4*x**2)**(sympy.S(7)/2)/448 + 9*(9 - 4*x**2)**(sympy.S(5)/2)/160 - 27*(9 - 4*x**2)**(sympy.S(3)/2)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_452():
    f = x**4*sqrt(9 - 4*x**2)
    F = x**5*sqrt(9 - 4*x**2)/6 - 3*x**3*sqrt(9 - 4*x**2)/32 - 81*x*sqrt(9 - 4*x**2)/256 + 729*asin(2*x/3)/512
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_453():
    f = x**3*sqrt(9 - 4*x**2)
    F = (9 - 4*x**2)**(sympy.S(5)/2)/80 - 3*(9 - 4*x**2)**(sympy.S(3)/2)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_454():
    f = x**2*sqrt(9 - 4*x**2)
    F = x**3*sqrt(9 - 4*x**2)/4 - 9*x*sqrt(9 - 4*x**2)/32 + 81*asin(2*x/3)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_455():
    f = x*sqrt(9 - 4*x**2)
    F = -(9 - 4*x**2)**(sympy.S(3)/2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_456():
    f = sqrt(9 - 4*x**2)
    F = x*sqrt(9 - 4*x**2)/2 + 9*asin(2*x/3)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_457():
    f = sqrt(9 - 4*x**2)/x
    F = sqrt(9 - 4*x**2) - 3*atanh(sqrt(9 - 4*x**2)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_458():
    f = sqrt(9 - 4*x**2)/x**2
    F = -2*asin(2*x/3) - sqrt(9 - 4*x**2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_459():
    f = sqrt(9 - 4*x**2)/x**3
    F = 2*atanh(sqrt(9 - 4*x**2)/3)/3 - sqrt(9 - 4*x**2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_460():
    f = sqrt(9 - 4*x**2)/x**4
    F = -(9 - 4*x**2)**(sympy.S(3)/2)/(27*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_461():
    f = sqrt(9 - 4*x**2)/x**5
    F = 2*atanh(sqrt(9 - 4*x**2)/3)/27 + sqrt(9 - 4*x**2)/(18*x**2) - sqrt(9 - 4*x**2)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_462():
    f = x**5*sqrt(4*x**2 - 9)
    F = (4*x**2 - 9)**(sympy.S(7)/2)/448 + 9*(4*x**2 - 9)**(sympy.S(5)/2)/160 + 27*(4*x**2 - 9)**(sympy.S(3)/2)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_463():
    f = x**4*sqrt(4*x**2 - 9)
    F = x**5*sqrt(4*x**2 - 9)/6 - 3*x**3*sqrt(4*x**2 - 9)/32 - 81*x*sqrt(4*x**2 - 9)/256 - 729*atanh(2*x/sqrt(4*x**2 - 9))/512
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_464():
    f = x**3*sqrt(4*x**2 - 9)
    F = (4*x**2 - 9)**(sympy.S(5)/2)/80 + 3*(4*x**2 - 9)**(sympy.S(3)/2)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_465():
    f = x**2*sqrt(4*x**2 - 9)
    F = x**3*sqrt(4*x**2 - 9)/4 - 9*x*sqrt(4*x**2 - 9)/32 - 81*atanh(2*x/sqrt(4*x**2 - 9))/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_466():
    f = x*sqrt(4*x**2 - 9)
    F = (4*x**2 - 9)**(sympy.S(3)/2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_467():
    f = sqrt(4*x**2 - 9)
    F = x*sqrt(4*x**2 - 9)/2 - 9*atanh(2*x/sqrt(4*x**2 - 9))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_468():
    f = sqrt(4*x**2 - 9)/x
    F = sqrt(4*x**2 - 9) - 3*atan(sqrt(4*x**2 - 9)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_469():
    f = sqrt(4*x**2 - 9)/x**2
    F = 2*atanh(2*x/sqrt(4*x**2 - 9)) - sqrt(4*x**2 - 9)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_470():
    f = sqrt(4*x**2 - 9)/x**3
    F = 2*atan(sqrt(4*x**2 - 9)/3)/3 - sqrt(4*x**2 - 9)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_471():
    f = sqrt(4*x**2 - 9)/x**4
    F = (4*x**2 - 9)**(sympy.S(3)/2)/(27*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_472():
    f = sqrt(4*x**2 - 9)/x**5
    F = 2*atan(sqrt(4*x**2 - 9)/3)/27 + sqrt(4*x**2 - 9)/(18*x**2) - sqrt(4*x**2 - 9)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_473():
    f = x**5*sqrt(-4*x**2 - 9)
    F = -(-4*x**2 - 9)**(sympy.S(7)/2)/448 - 9*(-4*x**2 - 9)**(sympy.S(5)/2)/160 - 27*(-4*x**2 - 9)**(sympy.S(3)/2)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_474():
    f = x**4*sqrt(-4*x**2 - 9)
    F = x**5*sqrt(-4*x**2 - 9)/6 + 3*x**3*sqrt(-4*x**2 - 9)/32 - 81*x*sqrt(-4*x**2 - 9)/256 - 729*atan(2*x/sqrt(-4*x**2 - 9))/512
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_475():
    f = x**3*sqrt(-4*x**2 - 9)
    F = (-4*x**2 - 9)**(sympy.S(5)/2)/80 + 3*(-4*x**2 - 9)**(sympy.S(3)/2)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_476():
    f = x**2*sqrt(-4*x**2 - 9)
    F = x**3*sqrt(-4*x**2 - 9)/4 + 9*x*sqrt(-4*x**2 - 9)/32 + 81*atan(2*x/sqrt(-4*x**2 - 9))/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_477():
    f = x*sqrt(-4*x**2 - 9)
    F = -(-4*x**2 - 9)**(sympy.S(3)/2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_478():
    f = sqrt(-4*x**2 - 9)
    F = x*sqrt(-4*x**2 - 9)/2 - 9*atan(2*x/sqrt(-4*x**2 - 9))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_479():
    f = sqrt(-4*x**2 - 9)/x
    F = sqrt(-4*x**2 - 9) - 3*atan(sqrt(-4*x**2 - 9)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_480():
    f = sqrt(-4*x**2 - 9)/x**2
    F = -2*atan(2*x/sqrt(-4*x**2 - 9)) - sqrt(-4*x**2 - 9)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_481():
    f = sqrt(-4*x**2 - 9)/x**3
    F = -2*atan(sqrt(-4*x**2 - 9)/3)/3 - sqrt(-4*x**2 - 9)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_482():
    f = sqrt(-4*x**2 - 9)/x**4
    F = (-4*x**2 - 9)**(sympy.S(3)/2)/(27*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_483():
    f = sqrt(-4*x**2 - 9)/x**5
    F = 2*atan(sqrt(-4*x**2 - 9)/3)/27 - sqrt(-4*x**2 - 9)/(18*x**2) - sqrt(-4*x**2 - 9)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_484():
    f = x**5/sqrt(a + b*x**2)
    F = a**2*sqrt(a + b*x**2)/b**3 - 2*a*(a + b*x**2)**(sympy.S(3)/2)/(3*b**3) + (a + b*x**2)**(sympy.S(5)/2)/(5*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_485():
    f = x**4/sqrt(a + b*x**2)
    F = 3*a**2*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(5)/2)) - 3*a*x*sqrt(a + b*x**2)/(8*b**2) + x**3*sqrt(a + b*x**2)/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_486():
    f = x**3/sqrt(a + b*x**2)
    F = -a*sqrt(a + b*x**2)/b**2 + (a + b*x**2)**(sympy.S(3)/2)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_487():
    f = x**2/sqrt(a + b*x**2)
    F = -a*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(3)/2)) + x*sqrt(a + b*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_488():
    f = x/sqrt(a + b*x**2)
    F = sqrt(a + b*x**2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_489():
    f = 1/sqrt(a + b*x**2)
    F = atanh(sqrt(b)*x/sqrt(a + b*x**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_490():
    f = 1/(x*sqrt(a + b*x**2))
    F = -atanh(sqrt(a + b*x**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_491():
    f = 1/(x**2*sqrt(a + b*x**2))
    F = -sqrt(a + b*x**2)/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_492():
    f = 1/(x**3*sqrt(a + b*x**2))
    F = -sqrt(a + b*x**2)/(2*a*x**2) + b*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_493():
    f = 1/(x**4*sqrt(a + b*x**2))
    F = -sqrt(a + b*x**2)/(3*a*x**3) + 2*b*sqrt(a + b*x**2)/(3*a**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_494():
    f = 1/(x**5*sqrt(a + b*x**2))
    F = -sqrt(a + b*x**2)/(4*a*x**4) + 3*b*sqrt(a + b*x**2)/(8*a**2*x**2) - 3*b**2*atanh(sqrt(a + b*x**2)/sqrt(a))/(8*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_495():
    f = x**5/(a + b*x**2)**(sympy.S(3)/2)
    F = -a**2/(b**3*sqrt(a + b*x**2)) - 2*a*sqrt(a + b*x**2)/b**3 + (a + b*x**2)**(sympy.S(3)/2)/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_496():
    f = x**4/(a + b*x**2)**(sympy.S(3)/2)
    F = -3*a*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(5)/2)) - x**3/(b*sqrt(a + b*x**2)) + 3*x*sqrt(a + b*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_497():
    f = x**3/(a + b*x**2)**(sympy.S(3)/2)
    F = a/(b**2*sqrt(a + b*x**2)) + sqrt(a + b*x**2)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_498():
    f = x**2/(a + b*x**2)**(sympy.S(3)/2)
    F = -x/(b*sqrt(a + b*x**2)) + atanh(sqrt(b)*x/sqrt(a + b*x**2))/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_499():
    f = x/(a + b*x**2)**(sympy.S(3)/2)
    F = -1/(b*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_500():
    f = (a + b*x**2)**(sympy.S(-3)/2)
    F = x/(a*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_501():
    f = 1/(x*(a + b*x**2)**(sympy.S(3)/2))
    F = 1/(a*sqrt(a + b*x**2)) - atanh(sqrt(a + b*x**2)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_502():
    f = 1/(x**2*(a + b*x**2)**(sympy.S(3)/2))
    F = -1/(a*x*sqrt(a + b*x**2)) - 2*b*x/(a**2*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_503():
    f = 1/(x**3*(a + b*x**2)**(sympy.S(3)/2))
    F = -1/(2*a*x**2*sqrt(a + b*x**2)) - 3*b/(2*a**2*sqrt(a + b*x**2)) + 3*b*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_504():
    f = 1/(x**4*(a + b*x**2)**(sympy.S(3)/2))
    F = -1/(3*a*x**3*sqrt(a + b*x**2)) + 4*b/(3*a**2*x*sqrt(a + b*x**2)) + 8*b**2*x/(3*a**3*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_505():
    f = x**6/(a + b*x**2)**(sympy.S(5)/2)
    F = -5*a*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(7)/2)) - x**5/(3*b*(a + b*x**2)**(sympy.S(3)/2)) - 5*x**3/(3*b**2*sqrt(a + b*x**2)) + 5*x*sqrt(a + b*x**2)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_506():
    f = x**5/(a + b*x**2)**(sympy.S(5)/2)
    F = -a**2/(3*b**3*(a + b*x**2)**(sympy.S(3)/2)) + 2*a/(b**3*sqrt(a + b*x**2)) + sqrt(a + b*x**2)/b**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_507():
    f = x**4/(a + b*x**2)**(sympy.S(5)/2)
    F = -x**3/(3*b*(a + b*x**2)**(sympy.S(3)/2)) - x/(b**2*sqrt(a + b*x**2)) + atanh(sqrt(b)*x/sqrt(a + b*x**2))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_508():
    f = x**3/(a + b*x**2)**(sympy.S(5)/2)
    F = a/(3*b**2*(a + b*x**2)**(sympy.S(3)/2)) - 1/(b**2*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_509():
    f = x**2/(a + b*x**2)**(sympy.S(5)/2)
    F = x**3/(3*a*(a + b*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_510():
    f = x/(a + b*x**2)**(sympy.S(5)/2)
    F = -1/(3*b*(a + b*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_511():
    f = (a + b*x**2)**(sympy.S(-5)/2)
    F = x/(3*a*(a + b*x**2)**(sympy.S(3)/2)) + 2*x/(3*a**2*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_512():
    f = 1/(x*(a + b*x**2)**(sympy.S(5)/2))
    F = 1/(3*a*(a + b*x**2)**(sympy.S(3)/2)) + 1/(a**2*sqrt(a + b*x**2)) - atanh(sqrt(a + b*x**2)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_513():
    f = 1/(x**2*(a + b*x**2)**(sympy.S(5)/2))
    F = -1/(a*x*(a + b*x**2)**(sympy.S(3)/2)) - 4*b*x/(3*a**2*(a + b*x**2)**(sympy.S(3)/2)) - 8*b*x/(3*a**3*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_514():
    f = 1/(x**3*(a + b*x**2)**(sympy.S(5)/2))
    F = -1/(2*a*x**2*(a + b*x**2)**(sympy.S(3)/2)) - 5*b/(6*a**2*(a + b*x**2)**(sympy.S(3)/2)) - 5*b/(2*a**3*sqrt(a + b*x**2)) + 5*b*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_515():
    f = 1/(x**4*(a + b*x**2)**(sympy.S(5)/2))
    F = -1/(3*a*x**3*(a + b*x**2)**(sympy.S(3)/2)) + 2*b/(a**2*x*(a + b*x**2)**(sympy.S(3)/2)) + 8*b**2*x/(3*a**3*(a + b*x**2)**(sympy.S(3)/2)) + 16*b**2*x/(3*a**4*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_516():
    f = x**10/(a + b*x**2)**(sympy.S(9)/2)
    F = -9*a*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(11)/2)) - x**9/(7*b*(a + b*x**2)**(sympy.S(7)/2)) - 9*x**7/(35*b**2*(a + b*x**2)**(sympy.S(5)/2)) - 3*x**5/(5*b**3*(a + b*x**2)**(sympy.S(3)/2)) - 3*x**3/(b**4*sqrt(a + b*x**2)) + 9*x*sqrt(a + b*x**2)/(2*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_517():
    f = x**9/(a + b*x**2)**(sympy.S(9)/2)
    F = -a**4/(7*b**5*(a + b*x**2)**(sympy.S(7)/2)) + 4*a**3/(5*b**5*(a + b*x**2)**(sympy.S(5)/2)) - 2*a**2/(b**5*(a + b*x**2)**(sympy.S(3)/2)) + 4*a/(b**5*sqrt(a + b*x**2)) + sqrt(a + b*x**2)/b**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_518():
    f = x**8/(a + b*x**2)**(sympy.S(9)/2)
    F = -x**7/(7*b*(a + b*x**2)**(sympy.S(7)/2)) - x**5/(5*b**2*(a + b*x**2)**(sympy.S(5)/2)) - x**3/(3*b**3*(a + b*x**2)**(sympy.S(3)/2)) - x/(b**4*sqrt(a + b*x**2)) + atanh(sqrt(b)*x/sqrt(a + b*x**2))/b**(sympy.S(9)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_519():
    f = x**7/(a + b*x**2)**(sympy.S(9)/2)
    F = a**3/(7*b**4*(a + b*x**2)**(sympy.S(7)/2)) - 3*a**2/(5*b**4*(a + b*x**2)**(sympy.S(5)/2)) + a/(b**4*(a + b*x**2)**(sympy.S(3)/2)) - 1/(b**4*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_520():
    f = x**6/(a + b*x**2)**(sympy.S(9)/2)
    F = x**7/(7*a*(a + b*x**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_521():
    f = x**5/(a + b*x**2)**(sympy.S(9)/2)
    F = -a**2/(7*b**3*(a + b*x**2)**(sympy.S(7)/2)) + 2*a/(5*b**3*(a + b*x**2)**(sympy.S(5)/2)) - 1/(3*b**3*(a + b*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_522():
    f = x**4/(a + b*x**2)**(sympy.S(9)/2)
    F = x**5/(5*a*(a + b*x**2)**(sympy.S(7)/2)) + 2*b*x**7/(35*a**2*(a + b*x**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_523():
    f = x**3/(a + b*x**2)**(sympy.S(9)/2)
    F = a/(7*b**2*(a + b*x**2)**(sympy.S(7)/2)) - 1/(5*b**2*(a + b*x**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_524():
    f = x**2/(a + b*x**2)**(sympy.S(9)/2)
    F = x**3/(3*a*(a + b*x**2)**(sympy.S(7)/2)) + 4*b*x**5/(15*a**2*(a + b*x**2)**(sympy.S(7)/2)) + 8*b**2*x**7/(105*a**3*(a + b*x**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_525():
    f = x/(a + b*x**2)**(sympy.S(9)/2)
    F = -1/(7*b*(a + b*x**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_526():
    f = (a + b*x**2)**(sympy.S(-9)/2)
    F = x/(7*a*(a + b*x**2)**(sympy.S(7)/2)) + 6*x/(35*a**2*(a + b*x**2)**(sympy.S(5)/2)) + 8*x/(35*a**3*(a + b*x**2)**(sympy.S(3)/2)) + 16*x/(35*a**4*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_527():
    f = 1/(x*(a + b*x**2)**(sympy.S(9)/2))
    F = 1/(7*a*(a + b*x**2)**(sympy.S(7)/2)) + 1/(5*a**2*(a + b*x**2)**(sympy.S(5)/2)) + 1/(3*a**3*(a + b*x**2)**(sympy.S(3)/2)) + 1/(a**4*sqrt(a + b*x**2)) - atanh(sqrt(a + b*x**2)/sqrt(a))/a**(sympy.S(9)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_528():
    f = 1/(x**2*(a + b*x**2)**(sympy.S(9)/2))
    F = -1/(a*x*(a + b*x**2)**(sympy.S(7)/2)) - 8*b*x/(7*a**2*(a + b*x**2)**(sympy.S(7)/2)) - 48*b*x/(35*a**3*(a + b*x**2)**(sympy.S(5)/2)) - 64*b*x/(35*a**4*(a + b*x**2)**(sympy.S(3)/2)) - 128*b*x/(35*a**5*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_529():
    f = 1/(x**3*(a + b*x**2)**(sympy.S(9)/2))
    F = -1/(2*a*x**2*(a + b*x**2)**(sympy.S(7)/2)) - 9*b/(14*a**2*(a + b*x**2)**(sympy.S(7)/2)) - 9*b/(10*a**3*(a + b*x**2)**(sympy.S(5)/2)) - 3*b/(2*a**4*(a + b*x**2)**(sympy.S(3)/2)) - 9*b/(2*a**5*sqrt(a + b*x**2)) + 9*b*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*a**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_530():
    f = 1/(x**4*(a + b*x**2)**(sympy.S(9)/2))
    F = -1/(3*a*x**3*(a + b*x**2)**(sympy.S(7)/2)) + 10*b/(3*a**2*x*(a + b*x**2)**(sympy.S(7)/2)) + 80*b**2*x/(21*a**3*(a + b*x**2)**(sympy.S(7)/2)) + 32*b**2*x/(7*a**4*(a + b*x**2)**(sympy.S(5)/2)) + 128*b**2*x/(21*a**5*(a + b*x**2)**(sympy.S(3)/2)) + 256*b**2*x/(21*a**6*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_531():
    f = x**5/sqrt(4*x**2 + 9)
    F = (4*x**2 + 9)**(sympy.S(5)/2)/320 - 3*(4*x**2 + 9)**(sympy.S(3)/2)/32 + 81*sqrt(4*x**2 + 9)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_532():
    f = x**4/sqrt(4*x**2 + 9)
    F = x**3*sqrt(4*x**2 + 9)/16 - 27*x*sqrt(4*x**2 + 9)/128 + 243*asinh(2*x/3)/256
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_533():
    f = x**3/sqrt(4*x**2 + 9)
    F = (4*x**2 + 9)**(sympy.S(3)/2)/48 - 9*sqrt(4*x**2 + 9)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_534():
    f = x**2/sqrt(4*x**2 + 9)
    F = x*sqrt(4*x**2 + 9)/8 - 9*asinh(2*x/3)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_535():
    f = x/sqrt(4*x**2 + 9)
    F = sqrt(4*x**2 + 9)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_536():
    f = 1/sqrt(4*x**2 + 9)
    F = asinh(2*x/3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_537():
    f = 1/(x*sqrt(4*x**2 + 9))
    F = -atanh(sqrt(4*x**2 + 9)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_538():
    f = 1/(x**2*sqrt(4*x**2 + 9))
    F = -sqrt(4*x**2 + 9)/(9*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_539():
    f = 1/(x**3*sqrt(4*x**2 + 9))
    F = 2*atanh(sqrt(4*x**2 + 9)/3)/27 - sqrt(4*x**2 + 9)/(18*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_540():
    f = 1/(x**4*sqrt(4*x**2 + 9))
    F = 8*sqrt(4*x**2 + 9)/(243*x) - sqrt(4*x**2 + 9)/(27*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_541():
    f = 1/(x**5*sqrt(4*x**2 + 9))
    F = -2*atanh(sqrt(4*x**2 + 9)/3)/81 + sqrt(4*x**2 + 9)/(54*x**2) - sqrt(4*x**2 + 9)/(36*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_542():
    f = x**5/sqrt(9 - 4*x**2)
    F = -(9 - 4*x**2)**(sympy.S(5)/2)/320 + 3*(9 - 4*x**2)**(sympy.S(3)/2)/32 - 81*sqrt(9 - 4*x**2)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_543():
    f = x**4/sqrt(9 - 4*x**2)
    F = -x**3*sqrt(9 - 4*x**2)/16 - 27*x*sqrt(9 - 4*x**2)/128 + 243*asin(2*x/3)/256
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_544():
    f = x**3/sqrt(9 - 4*x**2)
    F = (9 - 4*x**2)**(sympy.S(3)/2)/48 - 9*sqrt(9 - 4*x**2)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_545():
    f = x**2/sqrt(9 - 4*x**2)
    F = -x*sqrt(9 - 4*x**2)/8 + 9*asin(2*x/3)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_546():
    f = x/sqrt(9 - 4*x**2)
    F = -sqrt(9 - 4*x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_547():
    f = 1/sqrt(9 - 4*x**2)
    F = asin(2*x/3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_548():
    f = 1/(x*sqrt(9 - 4*x**2))
    F = -atanh(sqrt(9 - 4*x**2)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_549():
    f = 1/(x**2*sqrt(9 - 4*x**2))
    F = -sqrt(9 - 4*x**2)/(9*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_550():
    f = 1/(x**3*sqrt(9 - 4*x**2))
    F = -2*atanh(sqrt(9 - 4*x**2)/3)/27 - sqrt(9 - 4*x**2)/(18*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_551():
    f = 1/(x**4*sqrt(9 - 4*x**2))
    F = -8*sqrt(9 - 4*x**2)/(243*x) - sqrt(9 - 4*x**2)/(27*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_552():
    f = 1/(x**5*sqrt(9 - 4*x**2))
    F = -2*atanh(sqrt(9 - 4*x**2)/3)/81 - sqrt(9 - 4*x**2)/(54*x**2) - sqrt(9 - 4*x**2)/(36*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_553():
    f = x**5/sqrt(4*x**2 - 9)
    F = (4*x**2 - 9)**(sympy.S(5)/2)/320 + 3*(4*x**2 - 9)**(sympy.S(3)/2)/32 + 81*sqrt(4*x**2 - 9)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_554():
    f = x**4/sqrt(4*x**2 - 9)
    F = x**3*sqrt(4*x**2 - 9)/16 + 27*x*sqrt(4*x**2 - 9)/128 + 243*atanh(2*x/sqrt(4*x**2 - 9))/256
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_555():
    f = x**3/sqrt(4*x**2 - 9)
    F = (4*x**2 - 9)**(sympy.S(3)/2)/48 + 9*sqrt(4*x**2 - 9)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_556():
    f = x**2/sqrt(4*x**2 - 9)
    F = x*sqrt(4*x**2 - 9)/8 + 9*atanh(2*x/sqrt(4*x**2 - 9))/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_557():
    f = x/sqrt(4*x**2 - 9)
    F = sqrt(4*x**2 - 9)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_558():
    f = 1/sqrt(4*x**2 - 9)
    F = atanh(2*x/sqrt(4*x**2 - 9))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_559():
    f = 1/(x*sqrt(4*x**2 - 9))
    F = atan(sqrt(4*x**2 - 9)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_560():
    f = 1/(x**2*sqrt(4*x**2 - 9))
    F = sqrt(4*x**2 - 9)/(9*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_561():
    f = 1/(x**3*sqrt(4*x**2 - 9))
    F = 2*atan(sqrt(4*x**2 - 9)/3)/27 + sqrt(4*x**2 - 9)/(18*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_562():
    f = 1/(x**4*sqrt(4*x**2 - 9))
    F = 8*sqrt(4*x**2 - 9)/(243*x) + sqrt(4*x**2 - 9)/(27*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_563():
    f = 1/(x**5*sqrt(4*x**2 - 9))
    F = 2*atan(sqrt(4*x**2 - 9)/3)/81 + sqrt(4*x**2 - 9)/(54*x**2) + sqrt(4*x**2 - 9)/(36*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_564():
    f = x**5/sqrt(-4*x**2 - 9)
    F = -(-4*x**2 - 9)**(sympy.S(5)/2)/320 - 3*(-4*x**2 - 9)**(sympy.S(3)/2)/32 - 81*sqrt(-4*x**2 - 9)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_565():
    f = x**4/sqrt(-4*x**2 - 9)
    F = -x**3*sqrt(-4*x**2 - 9)/16 + 27*x*sqrt(-4*x**2 - 9)/128 + 243*atan(2*x/sqrt(-4*x**2 - 9))/256
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_566():
    f = x**3/sqrt(-4*x**2 - 9)
    F = (-4*x**2 - 9)**(sympy.S(3)/2)/48 + 9*sqrt(-4*x**2 - 9)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_567():
    f = x**2/sqrt(-4*x**2 - 9)
    F = -x*sqrt(-4*x**2 - 9)/8 - 9*atan(2*x/sqrt(-4*x**2 - 9))/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_568():
    f = x/sqrt(-4*x**2 - 9)
    F = -sqrt(-4*x**2 - 9)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_569():
    f = 1/sqrt(-4*x**2 - 9)
    F = atan(2*x/sqrt(-4*x**2 - 9))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_570():
    f = 1/(x*sqrt(-4*x**2 - 9))
    F = atan(sqrt(-4*x**2 - 9)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_571():
    f = 1/(x**2*sqrt(-4*x**2 - 9))
    F = sqrt(-4*x**2 - 9)/(9*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_572():
    f = 1/(x**3*sqrt(-4*x**2 - 9))
    F = -2*atan(sqrt(-4*x**2 - 9)/3)/27 + sqrt(-4*x**2 - 9)/(18*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_573():
    f = 1/(x**4*sqrt(-4*x**2 - 9))
    F = -8*sqrt(-4*x**2 - 9)/(243*x) + sqrt(-4*x**2 - 9)/(27*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_574():
    f = 1/(x**5*sqrt(-4*x**2 - 9))
    F = 2*atan(sqrt(-4*x**2 - 9)/3)/81 - sqrt(-4*x**2 - 9)/(54*x**2) + sqrt(-4*x**2 - 9)/(36*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_575():
    f = 1/sqrt(b*x**2 + 9)
    F = asinh(sqrt(b)*x/3)/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_576():
    f = 1/sqrt(-b*x**2 + 9)
    F = asin(sqrt(b)*x/3)/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_577():
    f = 1/sqrt(b*x**2 - 9)
    F = atanh(sqrt(b)*x/sqrt(b*x**2 - 9))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_578():
    f = 1/sqrt(-b*x**2 - 9)
    F = atan(sqrt(b)*x/sqrt(-b*x**2 - 9))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_579():
    f = 1/sqrt(b*x**2 + pi)
    F = asinh(sqrt(b)*x/sqrt(pi))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_580():
    f = 1/sqrt(-b*x**2 + pi)
    F = asin(sqrt(b)*x/sqrt(pi))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_581():
    f = 1/sqrt(b*x**2 - pi)
    F = atanh(sqrt(b)*x/sqrt(b*x**2 - pi))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_582():
    f = 1/sqrt(-b*x**2 - pi)
    F = atan(sqrt(b)*x/sqrt(-b*x**2 - pi))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_583():
    f = 1/sqrt(a + b*x**2)
    F = atanh(sqrt(b)*x/sqrt(a + b*x**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_584():
    f = 1/sqrt(a - b*x**2)
    F = atan(sqrt(b)*x/sqrt(a - b*x**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_585():
    f = 1/sqrt(-a + b*x**2)
    F = atanh(sqrt(b)*x/sqrt(-a + b*x**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_586():
    f = 1/sqrt(-a - b*x**2)
    F = atan(sqrt(b)*x/sqrt(-a - b*x**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_587():
    f = 1/sqrt(a**2 - x**2)
    F = atan(x/sqrt(a**2 - x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_588():
    f = (c*x)**(sympy.S(7)/2)*sqrt(a + b*x**2)
    F = 10*a**(sympy.S(11)/4)*c**(sympy.S(7)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(231*b**(sympy.S(9)/4)*sqrt(a + b*x**2)) - 20*a**2*c**3*sqrt(c*x)*sqrt(a + b*x**2)/(231*b**2) + 4*a*c*(c*x)**(sympy.S(5)/2)*sqrt(a + b*x**2)/(77*b) + 2*(c*x)**(sympy.S(9)/2)*sqrt(a + b*x**2)/(11*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_589():
    f = (c*x)**(sympy.S(5)/2)*sqrt(a + b*x**2)
    F = 4*a**(sympy.S(9)/4)*c**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(15*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) - 2*a**(sympy.S(9)/4)*c**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(15*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) - 4*a**2*c**2*sqrt(c*x)*sqrt(a + b*x**2)/(15*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x)) + 4*a*c*(c*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)/(45*b) + 2*(c*x)**(sympy.S(7)/2)*sqrt(a + b*x**2)/(9*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_590():
    f = (c*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)
    F = -2*a**(sympy.S(7)/4)*c**(sympy.S(3)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(21*b**(sympy.S(5)/4)*sqrt(a + b*x**2)) + 4*a*c*sqrt(c*x)*sqrt(a + b*x**2)/(21*b) + 2*(c*x)**(sympy.S(5)/2)*sqrt(a + b*x**2)/(7*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_591():
    f = sqrt(c*x)*sqrt(a + b*x**2)
    F = -4*a**(sympy.S(5)/4)*sqrt(c)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(a + b*x**2)) + 2*a**(sympy.S(5)/4)*sqrt(c)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(a + b*x**2)) + 4*a*sqrt(c*x)*sqrt(a + b*x**2)/(5*sqrt(b)*(sqrt(a) + sqrt(b)*x)) + 2*(c*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)/(5*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_592():
    f = sqrt(a + b*x**2)/sqrt(c*x)
    F = 2*a**(sympy.S(3)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(3*b**(sympy.S(1)/4)*sqrt(c)*sqrt(a + b*x**2)) + 2*sqrt(c*x)*sqrt(a + b*x**2)/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_593():
    f = sqrt(a + b*x**2)/(c*x)**(sympy.S(3)/2)
    F = -4*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(c**(sympy.S(3)/2)*sqrt(a + b*x**2)) + 2*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(c**(sympy.S(3)/2)*sqrt(a + b*x**2)) + 4*sqrt(b)*sqrt(c*x)*sqrt(a + b*x**2)/(c**2*(sqrt(a) + sqrt(b)*x)) - 2*sqrt(a + b*x**2)/(c*sqrt(c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_594():
    f = sqrt(a + b*x**2)/(c*x)**(sympy.S(5)/2)
    F = -2*sqrt(a + b*x**2)/(3*c*(c*x)**(sympy.S(3)/2)) + 2*b**(sympy.S(3)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(3*a**(sympy.S(1)/4)*c**(sympy.S(5)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_595():
    f = sqrt(a + b*x**2)/(c*x)**(sympy.S(7)/2)
    F = -2*sqrt(a + b*x**2)/(5*c*(c*x)**(sympy.S(5)/2)) + 4*b**(sympy.S(3)/2)*sqrt(c*x)*sqrt(a + b*x**2)/(5*a*c**4*(sqrt(a) + sqrt(b)*x)) - 4*b*sqrt(a + b*x**2)/(5*a*c**3*sqrt(c*x)) - 4*b**(sympy.S(5)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(5*a**(sympy.S(3)/4)*c**(sympy.S(7)/2)*sqrt(a + b*x**2)) + 2*b**(sympy.S(5)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(5*a**(sympy.S(3)/4)*c**(sympy.S(7)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_596():
    f = (c*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(3)/2)
    F = 4*a**(sympy.S(15)/4)*c**(sympy.S(7)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(231*b**(sympy.S(9)/4)*sqrt(a + b*x**2)) - 8*a**3*c**3*sqrt(c*x)*sqrt(a + b*x**2)/(231*b**2) + 8*a**2*c*(c*x)**(sympy.S(5)/2)*sqrt(a + b*x**2)/(385*b) + 4*a*(c*x)**(sympy.S(9)/2)*sqrt(a + b*x**2)/(55*c) + 2*(c*x)**(sympy.S(9)/2)*(a + b*x**2)**(sympy.S(3)/2)/(15*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_597():
    f = (c*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/2)
    F = 8*a**(sympy.S(13)/4)*c**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(65*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) - 4*a**(sympy.S(13)/4)*c**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(65*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) - 8*a**3*c**2*sqrt(c*x)*sqrt(a + b*x**2)/(65*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x)) + 8*a**2*c*(c*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)/(195*b) + 4*a*(c*x)**(sympy.S(7)/2)*sqrt(a + b*x**2)/(39*c) + 2*(c*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(3)/2)/(13*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_598():
    f = (c*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/2)
    F = -4*a**(sympy.S(11)/4)*c**(sympy.S(3)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(77*b**(sympy.S(5)/4)*sqrt(a + b*x**2)) + 8*a**2*c*sqrt(c*x)*sqrt(a + b*x**2)/(77*b) + 12*a*(c*x)**(sympy.S(5)/2)*sqrt(a + b*x**2)/(77*c) + 2*(c*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/2)/(11*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_599():
    f = sqrt(c*x)*(a + b*x**2)**(sympy.S(3)/2)
    F = -8*a**(sympy.S(9)/4)*sqrt(c)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(15*b**(sympy.S(3)/4)*sqrt(a + b*x**2)) + 4*a**(sympy.S(9)/4)*sqrt(c)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(15*b**(sympy.S(3)/4)*sqrt(a + b*x**2)) + 8*a**2*sqrt(c*x)*sqrt(a + b*x**2)/(15*sqrt(b)*(sqrt(a) + sqrt(b)*x)) + 4*a*(c*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)/(15*c) + 2*(c*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/2)/(9*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_600():
    f = (a + b*x**2)**(sympy.S(3)/2)/sqrt(c*x)
    F = 4*a**(sympy.S(7)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(7*b**(sympy.S(1)/4)*sqrt(c)*sqrt(a + b*x**2)) + 4*a*sqrt(c*x)*sqrt(a + b*x**2)/(7*c) + 2*sqrt(c*x)*(a + b*x**2)**(sympy.S(3)/2)/(7*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_601():
    f = (a + b*x**2)**(sympy.S(3)/2)/(c*x)**(sympy.S(3)/2)
    F = -24*a**(sympy.S(5)/4)*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(5*c**(sympy.S(3)/2)*sqrt(a + b*x**2)) + 12*a**(sympy.S(5)/4)*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(5*c**(sympy.S(3)/2)*sqrt(a + b*x**2)) + 24*a*sqrt(b)*sqrt(c*x)*sqrt(a + b*x**2)/(5*c**2*(sqrt(a) + sqrt(b)*x)) + 12*b*(c*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)/(5*c**3) - 2*(a + b*x**2)**(sympy.S(3)/2)/(c*sqrt(c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_602():
    f = (a + b*x**2)**(sympy.S(3)/2)/(c*x)**(sympy.S(5)/2)
    F = 4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(3*c**(sympy.S(5)/2)*sqrt(a + b*x**2)) + 4*b*sqrt(c*x)*sqrt(a + b*x**2)/(3*c**3) - 2*(a + b*x**2)**(sympy.S(3)/2)/(3*c*(c*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_603():
    f = (a + b*x**2)**(sympy.S(3)/2)/(c*x)**(sympy.S(7)/2)
    F = -24*a**(sympy.S(1)/4)*b**(sympy.S(5)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(5*c**(sympy.S(7)/2)*sqrt(a + b*x**2)) + 12*a**(sympy.S(1)/4)*b**(sympy.S(5)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(5*c**(sympy.S(7)/2)*sqrt(a + b*x**2)) + 24*b**(sympy.S(3)/2)*sqrt(c*x)*sqrt(a + b*x**2)/(5*c**4*(sqrt(a) + sqrt(b)*x)) - 12*b*sqrt(a + b*x**2)/(5*c**3*sqrt(c*x)) - 2*(a + b*x**2)**(sympy.S(3)/2)/(5*c*(c*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_604():
    f = (a + b*x**2)**(sympy.S(3)/2)/(c*x)**(sympy.S(9)/2)
    F = -4*b*sqrt(a + b*x**2)/(7*c**3*(c*x)**(sympy.S(3)/2)) - 2*(a + b*x**2)**(sympy.S(3)/2)/(7*c*(c*x)**(sympy.S(7)/2)) + 4*b**(sympy.S(7)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(7*a**(sympy.S(1)/4)*c**(sympy.S(9)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_605():
    f = (a + b*x**2)**(sympy.S(3)/2)/(c*x)**(sympy.S(11)/2)
    F = -4*b*sqrt(a + b*x**2)/(15*c**3*(c*x)**(sympy.S(5)/2)) - 2*(a + b*x**2)**(sympy.S(3)/2)/(9*c*(c*x)**(sympy.S(9)/2)) + 8*b**(sympy.S(5)/2)*sqrt(c*x)*sqrt(a + b*x**2)/(15*a*c**6*(sqrt(a) + sqrt(b)*x)) - 8*b**2*sqrt(a + b*x**2)/(15*a*c**5*sqrt(c*x)) - 8*b**(sympy.S(9)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(15*a**(sympy.S(3)/4)*c**(sympy.S(11)/2)*sqrt(a + b*x**2)) + 4*b**(sympy.S(9)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(15*a**(sympy.S(3)/4)*c**(sympy.S(11)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_606():
    f = (c*x)**(sympy.S(5)/2)*sqrt(-2*a*x**2 + 3*a)
    F = -3*6**(sympy.S(1)/4)*a*c**2*sqrt(c*x)*sqrt(3 - 2*x**2)*elliptic_e(asin(sqrt(6)*sqrt(-sqrt(6)*x + 3)/6), 2)/(5*sqrt(x)*sqrt(-2*a*x**2 + 3*a)) - 2*c*(c*x)**(sympy.S(3)/2)*sqrt(-2*a*x**2 + 3*a)/15 + 2*(c*x)**(sympy.S(7)/2)*sqrt(-2*a*x**2 + 3*a)/(9*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_607():
    f = (c*x)**(sympy.S(3)/2)*sqrt(-2*a*x**2 + 3*a)
    F = 6**(sympy.S(3)/4)*a*c**(sympy.S(3)/2)*sqrt(3 - 2*x**2)*elliptic_f(asin(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*sqrt(c*x)/(3*sqrt(c))), -1)/(7*sqrt(a*(3 - 2*x**2))) - 2*c*sqrt(c*x)*sqrt(-2*a*x**2 + 3*a)/7 + 2*(c*x)**(sympy.S(5)/2)*sqrt(-2*a*x**2 + 3*a)/(7*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_608():
    f = sqrt(c*x)*sqrt(-2*a*x**2 + 3*a)
    F = -6*6**(sympy.S(1)/4)*a*sqrt(c*x)*sqrt(3 - 2*x**2)*elliptic_e(asin(sqrt(6)*sqrt(-sqrt(6)*x + 3)/6), 2)/(5*sqrt(x)*sqrt(-2*a*x**2 + 3*a)) + 2*(c*x)**(sympy.S(3)/2)*sqrt(-2*a*x**2 + 3*a)/(5*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_609():
    f = sqrt(-2*a*x**2 + 3*a)/sqrt(c*x)
    F = 2*6**(sympy.S(3)/4)*a*sqrt(3 - 2*x**2)*elliptic_f(asin(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*sqrt(c*x)/(3*sqrt(c))), -1)/(3*sqrt(c)*sqrt(a*(3 - 2*x**2))) + 2*sqrt(c*x)*sqrt(-2*a*x**2 + 3*a)/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_610():
    f = sqrt(-2*a*x**2 + 3*a)/(c*x)**(sympy.S(3)/2)
    F = 4*6**(sympy.S(1)/4)*a*sqrt(c*x)*sqrt(3 - 2*x**2)*elliptic_e(asin(sqrt(6)*sqrt(-sqrt(6)*x + 3)/6), 2)/(c**2*sqrt(x)*sqrt(-2*a*x**2 + 3*a)) - 2*sqrt(-2*a*x**2 + 3*a)/(c*sqrt(c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_611():
    f = sqrt(-2*a*x**2 + 3*a)/(c*x)**(sympy.S(5)/2)
    F = -4*6**(sympy.S(3)/4)*a*sqrt(3 - 2*x**2)*elliptic_f(asin(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*sqrt(c*x)/(3*sqrt(c))), -1)/(9*c**(sympy.S(5)/2)*sqrt(a*(3 - 2*x**2))) - 2*sqrt(-2*a*x**2 + 3*a)/(3*c*(c*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_612():
    f = (c*x)**(sympy.S(7)/2)/sqrt(a + b*x**2)
    F = 5*a**(sympy.S(7)/4)*c**(sympy.S(7)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(21*b**(sympy.S(9)/4)*sqrt(a + b*x**2)) - 10*a*c**3*sqrt(c*x)*sqrt(a + b*x**2)/(21*b**2) + 2*c*(c*x)**(sympy.S(5)/2)*sqrt(a + b*x**2)/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_613():
    f = (c*x)**(sympy.S(5)/2)/sqrt(a + b*x**2)
    F = 6*a**(sympy.S(5)/4)*c**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(5*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) - 3*a**(sympy.S(5)/4)*c**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(5*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) - 6*a*c**2*sqrt(c*x)*sqrt(a + b*x**2)/(5*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x)) + 2*c*(c*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_614():
    f = (c*x)**(sympy.S(3)/2)/sqrt(a + b*x**2)
    F = -a**(sympy.S(3)/4)*c**(sympy.S(3)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(3*b**(sympy.S(5)/4)*sqrt(a + b*x**2)) + 2*c*sqrt(c*x)*sqrt(a + b*x**2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_615():
    f = sqrt(c*x)/sqrt(a + b*x**2)
    F = -2*a**(sympy.S(1)/4)*sqrt(c)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(b**(sympy.S(3)/4)*sqrt(a + b*x**2)) + a**(sympy.S(1)/4)*sqrt(c)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(b**(sympy.S(3)/4)*sqrt(a + b*x**2)) + 2*sqrt(c*x)*sqrt(a + b*x**2)/(sqrt(b)*(sqrt(a) + sqrt(b)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_616():
    f = 1/(sqrt(c*x)*sqrt(a + b*x**2))
    F = sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(c)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_617():
    f = 1/((c*x)**(sympy.S(3)/2)*sqrt(a + b*x**2))
    F = 2*sqrt(b)*sqrt(c*x)*sqrt(a + b*x**2)/(a*c**2*(sqrt(a) + sqrt(b)*x)) - 2*sqrt(a + b*x**2)/(a*c*sqrt(c*x)) - 2*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(a**(sympy.S(3)/4)*c**(sympy.S(3)/2)*sqrt(a + b*x**2)) + b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(a**(sympy.S(3)/4)*c**(sympy.S(3)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_618():
    f = 1/((c*x)**(sympy.S(5)/2)*sqrt(a + b*x**2))
    F = -2*sqrt(a + b*x**2)/(3*a*c*(c*x)**(sympy.S(3)/2)) - b**(sympy.S(3)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(3*a**(sympy.S(5)/4)*c**(sympy.S(5)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_619():
    f = 1/((c*x)**(sympy.S(7)/2)*sqrt(a + b*x**2))
    F = -2*sqrt(a + b*x**2)/(5*a*c*(c*x)**(sympy.S(5)/2)) - 6*b**(sympy.S(3)/2)*sqrt(c*x)*sqrt(a + b*x**2)/(5*a**2*c**4*(sqrt(a) + sqrt(b)*x)) + 6*b*sqrt(a + b*x**2)/(5*a**2*c**3*sqrt(c*x)) + 6*b**(sympy.S(5)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(5*a**(sympy.S(7)/4)*c**(sympy.S(7)/2)*sqrt(a + b*x**2)) - 3*b**(sympy.S(5)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(5*a**(sympy.S(7)/4)*c**(sympy.S(7)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_620():
    f = (c*x)**(sympy.S(7)/2)/(a + b*x**2)**(sympy.S(3)/2)
    F = -5*a**(sympy.S(3)/4)*c**(sympy.S(7)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(6*b**(sympy.S(9)/4)*sqrt(a + b*x**2)) - c*(c*x)**(sympy.S(5)/2)/(b*sqrt(a + b*x**2)) + 5*c**3*sqrt(c*x)*sqrt(a + b*x**2)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_621():
    f = (c*x)**(sympy.S(5)/2)/(a + b*x**2)**(sympy.S(3)/2)
    F = -3*a**(sympy.S(1)/4)*c**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(b**(sympy.S(7)/4)*sqrt(a + b*x**2)) + 3*a**(sympy.S(1)/4)*c**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(2*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) - c*(c*x)**(sympy.S(3)/2)/(b*sqrt(a + b*x**2)) + 3*c**2*sqrt(c*x)*sqrt(a + b*x**2)/(b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_622():
    f = (c*x)**(sympy.S(3)/2)/(a + b*x**2)**(sympy.S(3)/2)
    F = -c*sqrt(c*x)/(b*sqrt(a + b*x**2)) + c**(sympy.S(3)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(2*a**(sympy.S(1)/4)*b**(sympy.S(5)/4)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_623():
    f = sqrt(c*x)/(a + b*x**2)**(sympy.S(3)/2)
    F = (c*x)**(sympy.S(3)/2)/(a*c*sqrt(a + b*x**2)) - sqrt(c*x)*sqrt(a + b*x**2)/(a*sqrt(b)*(sqrt(a) + sqrt(b)*x)) + sqrt(c)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*sqrt(a + b*x**2)) - sqrt(c)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_624():
    f = 1/(sqrt(c*x)*(a + b*x**2)**(sympy.S(3)/2))
    F = sqrt(c*x)/(a*c*sqrt(a + b*x**2)) + sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(2*a**(sympy.S(5)/4)*b**(sympy.S(1)/4)*sqrt(c)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_625():
    f = 1/((c*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/2))
    F = 1/(a*c*sqrt(c*x)*sqrt(a + b*x**2)) + 3*sqrt(b)*sqrt(c*x)*sqrt(a + b*x**2)/(a**2*c**2*(sqrt(a) + sqrt(b)*x)) - 3*sqrt(a + b*x**2)/(a**2*c*sqrt(c*x)) - 3*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(a**(sympy.S(7)/4)*c**(sympy.S(3)/2)*sqrt(a + b*x**2)) + 3*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(2*a**(sympy.S(7)/4)*c**(sympy.S(3)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_626():
    f = 1/((c*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/2))
    F = 1/(a*c*(c*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)) - 5*sqrt(a + b*x**2)/(3*a**2*c*(c*x)**(sympy.S(3)/2)) - 5*b**(sympy.S(3)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(6*a**(sympy.S(9)/4)*c**(sympy.S(5)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_627():
    f = 1/((c*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(3)/2))
    F = 1/(a*c*(c*x)**(sympy.S(5)/2)*sqrt(a + b*x**2)) - 7*sqrt(a + b*x**2)/(5*a**2*c*(c*x)**(sympy.S(5)/2)) - 21*b**(sympy.S(3)/2)*sqrt(c*x)*sqrt(a + b*x**2)/(5*a**3*c**4*(sqrt(a) + sqrt(b)*x)) + 21*b*sqrt(a + b*x**2)/(5*a**3*c**3*sqrt(c*x)) + 21*b**(sympy.S(5)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(5*a**(sympy.S(11)/4)*c**(sympy.S(7)/2)*sqrt(a + b*x**2)) - 21*b**(sympy.S(5)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(10*a**(sympy.S(11)/4)*c**(sympy.S(7)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_628():
    f = (c*x)**(sympy.S(7)/2)/(a + b*x**2)**(sympy.S(5)/2)
    F = -c*(c*x)**(sympy.S(5)/2)/(3*b*(a + b*x**2)**(sympy.S(3)/2)) - 5*c**3*sqrt(c*x)/(6*b**2*sqrt(a + b*x**2)) + 5*c**(sympy.S(7)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(12*a**(sympy.S(1)/4)*b**(sympy.S(9)/4)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_629():
    f = (c*x)**(sympy.S(5)/2)/(a + b*x**2)**(sympy.S(5)/2)
    F = -c*(c*x)**(sympy.S(3)/2)/(3*b*(a + b*x**2)**(sympy.S(3)/2)) + c*(c*x)**(sympy.S(3)/2)/(2*a*b*sqrt(a + b*x**2)) - c**2*sqrt(c*x)*sqrt(a + b*x**2)/(2*a*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x)) + c**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(2*a**(sympy.S(3)/4)*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) - c**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(4*a**(sympy.S(3)/4)*b**(sympy.S(7)/4)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_630():
    f = (c*x)**(sympy.S(3)/2)/(a + b*x**2)**(sympy.S(5)/2)
    F = -c*sqrt(c*x)/(3*b*(a + b*x**2)**(sympy.S(3)/2)) + c*sqrt(c*x)/(6*a*b*sqrt(a + b*x**2)) + c**(sympy.S(3)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(12*a**(sympy.S(5)/4)*b**(sympy.S(5)/4)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_631():
    f = sqrt(c*x)/(a + b*x**2)**(sympy.S(5)/2)
    F = (c*x)**(sympy.S(3)/2)/(3*a*c*(a + b*x**2)**(sympy.S(3)/2)) + (c*x)**(sympy.S(3)/2)/(2*a**2*c*sqrt(a + b*x**2)) - sqrt(c*x)*sqrt(a + b*x**2)/(2*a**2*sqrt(b)*(sqrt(a) + sqrt(b)*x)) + sqrt(c)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(2*a**(sympy.S(7)/4)*b**(sympy.S(3)/4)*sqrt(a + b*x**2)) - sqrt(c)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(4*a**(sympy.S(7)/4)*b**(sympy.S(3)/4)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_632():
    f = 1/(sqrt(c*x)*(a + b*x**2)**(sympy.S(5)/2))
    F = sqrt(c*x)/(3*a*c*(a + b*x**2)**(sympy.S(3)/2)) + 5*sqrt(c*x)/(6*a**2*c*sqrt(a + b*x**2)) + sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(5*sqrt(a) + 5*sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(12*a**(sympy.S(9)/4)*b**(sympy.S(1)/4)*sqrt(c)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_633():
    f = 1/((c*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(5)/2))
    F = 1/(3*a*c*sqrt(c*x)*(a + b*x**2)**(sympy.S(3)/2)) + 7/(6*a**2*c*sqrt(c*x)*sqrt(a + b*x**2)) + 7*sqrt(b)*sqrt(c*x)*sqrt(a + b*x**2)/(2*a**3*c**2*(sqrt(a) + sqrt(b)*x)) - 7*sqrt(a + b*x**2)/(2*a**3*c*sqrt(c*x)) - 7*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(2*a**(sympy.S(11)/4)*c**(sympy.S(3)/2)*sqrt(a + b*x**2)) + 7*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(4*a**(sympy.S(11)/4)*c**(sympy.S(3)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_634():
    f = 1/((c*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(5)/2))
    F = 1/(3*a*c*(c*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/2)) + 3/(2*a**2*c*(c*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)) - 5*sqrt(a + b*x**2)/(2*a**3*c*(c*x)**(sympy.S(3)/2)) - 5*b**(sympy.S(3)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(4*a**(sympy.S(13)/4)*c**(sympy.S(5)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_635():
    f = 1/((c*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(5)/2))
    F = 1/(3*a*c*(c*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/2)) + 11/(6*a**2*c*(c*x)**(sympy.S(5)/2)*sqrt(a + b*x**2)) - 77*sqrt(a + b*x**2)/(30*a**3*c*(c*x)**(sympy.S(5)/2)) - 77*b**(sympy.S(3)/2)*sqrt(c*x)*sqrt(a + b*x**2)/(10*a**4*c**4*(sqrt(a) + sqrt(b)*x)) + 77*b*sqrt(a + b*x**2)/(10*a**4*c**3*sqrt(c*x)) + 77*b**(sympy.S(5)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(10*a**(sympy.S(15)/4)*c**(sympy.S(7)/2)*sqrt(a + b*x**2)) - 77*b**(sympy.S(5)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(a**(sympy.S(1)/4)*sqrt(c))), sympy.S.Half)/(20*a**(sympy.S(15)/4)*c**(sympy.S(7)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_636():
    f = (c*x)**(sympy.S(5)/2)/sqrt(-2*a*x**2 + 3*a)
    F = -9*6**(sympy.S(1)/4)*c**2*sqrt(c*x)*sqrt(3 - 2*x**2)*elliptic_e(asin(sqrt(6)*sqrt(-sqrt(6)*x + 3)/6), 2)/(10*sqrt(x)*sqrt(-2*a*x**2 + 3*a)) - c*(c*x)**(sympy.S(3)/2)*sqrt(-2*a*x**2 + 3*a)/(5*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_637():
    f = (c*x)**(sympy.S(3)/2)/sqrt(-2*a*x**2 + 3*a)
    F = 6**(sympy.S(3)/4)*c**(sympy.S(3)/2)*sqrt(3 - 2*x**2)*elliptic_f(asin(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*sqrt(c*x)/(3*sqrt(c))), -1)/(6*sqrt(a*(3 - 2*x**2))) - c*sqrt(c*x)*sqrt(-2*a*x**2 + 3*a)/(3*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_638():
    f = sqrt(c*x)/sqrt(-2*a*x**2 + 3*a)
    F = -6**(sympy.S(1)/4)*sqrt(c*x)*sqrt(3 - 2*x**2)*elliptic_e(asin(sqrt(6)*sqrt(-sqrt(6)*x + 3)/6), 2)/(sqrt(x)*sqrt(-2*a*x**2 + 3*a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_639():
    f = 1/(sqrt(c*x)*sqrt(-2*a*x**2 + 3*a))
    F = 6**(sympy.S(3)/4)*sqrt(3 - 2*x**2)*elliptic_f(asin(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*sqrt(c*x)/(3*sqrt(c))), -1)/(3*sqrt(c)*sqrt(a*(3 - 2*x**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_640():
    f = 1/((c*x)**(sympy.S(3)/2)*sqrt(-2*a*x**2 + 3*a))
    F = 2*6**(sympy.S(1)/4)*sqrt(c*x)*sqrt(3 - 2*x**2)*elliptic_e(asin(sqrt(6)*sqrt(-sqrt(6)*x + 3)/6), 2)/(3*c**2*sqrt(x)*sqrt(-2*a*x**2 + 3*a)) - 2*sqrt(-2*a*x**2 + 3*a)/(3*a*c*sqrt(c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_641():
    f = 1/((c*x)**(sympy.S(5)/2)*sqrt(-2*a*x**2 + 3*a))
    F = 2*6**(sympy.S(3)/4)*sqrt(3 - 2*x**2)*elliptic_f(asin(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*sqrt(c*x)/(3*sqrt(c))), -1)/(27*c**(sympy.S(5)/2)*sqrt(a*(3 - 2*x**2))) - 2*sqrt(-2*a*x**2 + 3*a)/(9*a*c*(c*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_642():
    f = (c*x)**(sympy.S(5)/2)/(-2*a*x**2 + 3*a)**(sympy.S(3)/2)
    F = 3*6**(sympy.S(1)/4)*c**2*sqrt(c*x)*sqrt(3 - 2*x**2)*elliptic_e(asin(sqrt(6)*sqrt(-sqrt(6)*x + 3)/6), 2)/(4*a*sqrt(x)*sqrt(-2*a*x**2 + 3*a)) + c*(c*x)**(sympy.S(3)/2)/(2*a*sqrt(-2*a*x**2 + 3*a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_643():
    f = (c*x)**(sympy.S(3)/2)/(-2*a*x**2 + 3*a)**(sympy.S(3)/2)
    F = -6**(sympy.S(3)/4)*c**(sympy.S(3)/2)*sqrt(3 - 2*x**2)*elliptic_f(asin(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*sqrt(c*x)/(3*sqrt(c))), -1)/(12*a*sqrt(a*(3 - 2*x**2))) + c*sqrt(c*x)/(2*a*sqrt(-2*a*x**2 + 3*a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_644():
    f = sqrt(c*x)/(-2*a*x**2 + 3*a)**(sympy.S(3)/2)
    F = 6**(sympy.S(1)/4)*sqrt(c*x)*sqrt(3 - 2*x**2)*elliptic_e(asin(sqrt(6)*sqrt(-sqrt(6)*x + 3)/6), 2)/(6*a*sqrt(x)*sqrt(-2*a*x**2 + 3*a)) + (c*x)**(sympy.S(3)/2)/(3*a*c*sqrt(-2*a*x**2 + 3*a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_645():
    f = 1/(sqrt(c*x)*(-2*a*x**2 + 3*a)**(sympy.S(3)/2))
    F = sqrt(c*x)/(3*a*c*sqrt(-2*a*x**2 + 3*a)) + 6**(sympy.S(3)/4)*sqrt(3 - 2*x**2)*elliptic_f(asin(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*sqrt(c*x)/(3*sqrt(c))), -1)/(18*a*sqrt(c)*sqrt(a*(3 - 2*x**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_646():
    f = 1/((c*x)**(sympy.S(3)/2)*(-2*a*x**2 + 3*a)**(sympy.S(3)/2))
    F = 1/(3*a*c*sqrt(c*x)*sqrt(-2*a*x**2 + 3*a)) + 6**(sympy.S(1)/4)*sqrt(c*x)*sqrt(3 - 2*x**2)*elliptic_e(asin(sqrt(6)*sqrt(-sqrt(6)*x + 3)/6), 2)/(3*a*c**2*sqrt(x)*sqrt(-2*a*x**2 + 3*a)) - sqrt(-2*a*x**2 + 3*a)/(3*a**2*c*sqrt(c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_647():
    f = 1/((c*x)**(sympy.S(5)/2)*(-2*a*x**2 + 3*a)**(sympy.S(3)/2))
    F = 1/(3*a*c*(c*x)**(sympy.S(3)/2)*sqrt(-2*a*x**2 + 3*a)) + 5*6**(sympy.S(3)/4)*sqrt(3 - 2*x**2)*elliptic_f(asin(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*sqrt(c*x)/(3*sqrt(c))), -1)/(81*a*c**(sympy.S(5)/2)*sqrt(a*(3 - 2*x**2))) - 5*sqrt(-2*a*x**2 + 3*a)/(27*a**2*c*(c*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_648():
    f = 1/(sqrt(x)*sqrt(-a**2*x**2 + 1))
    F = 2*elliptic_f(asin(sqrt(a)*sqrt(x)), -1)/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_649():
    f = 1/(sqrt(x)*sqrt(a*x**2 + 1))
    F = sqrt((a*x**2 + 1)/(sqrt(a)*x + 1)**2)*(sqrt(a)*x + 1)*elliptic_f(2*atan(a**(sympy.S(1)/4)*sqrt(x)), sympy.S.Half)/(a**(sympy.S(1)/4)*sqrt(a*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_650():
    f = x**(m + 1)*(a*(m + 2) + b*x**2*(m + 3))/sqrt(a + b*x**2)
    F = x**(m + 2)*sqrt(a + b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_651():
    f = a*x**(m + 1)*(m + 2)/sqrt(a + b*x**2) + b*x**(m + 3)*(m + 3)/sqrt(a + b*x**2)
    F = x**(m + 2)*sqrt(a + b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_652():
    f = x**(m - 1)*(a*m + b*x**2*(m - 1))/(a + b*x**2)**(sympy.S(3)/2)
    F = x**m/sqrt(a + b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_653():
    f = x**7*(a + b*x**2)**(sympy.S(1)/3)
    F = -3*a**3*(a + b*x**2)**(sympy.S(4)/3)/(8*b**4) + 9*a**2*(a + b*x**2)**(sympy.S(7)/3)/(14*b**4) - 9*a*(a + b*x**2)**(sympy.S(10)/3)/(20*b**4) + 3*(a + b*x**2)**(sympy.S(13)/3)/(26*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_654():
    f = x**5*(a + b*x**2)**(sympy.S(1)/3)
    F = 3*a**2*(a + b*x**2)**(sympy.S(4)/3)/(8*b**3) - 3*a*(a + b*x**2)**(sympy.S(7)/3)/(7*b**3) + 3*(a + b*x**2)**(sympy.S(10)/3)/(20*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_655():
    f = x**3*(a + b*x**2)**(sympy.S(1)/3)
    F = -3*a*(a + b*x**2)**(sympy.S(4)/3)/(8*b**2) + 3*(a + b*x**2)**(sympy.S(7)/3)/(14*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_656():
    f = x*(a + b*x**2)**(sympy.S(1)/3)
    F = 3*(a + b*x**2)**(sympy.S(4)/3)/(8*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_657():
    f = (a + b*x**2)**(sympy.S(1)/3)/x
    F = -a**(sympy.S(1)/3)*log(x)/2 + 3*a**(sympy.S(1)/3)*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/4 - sqrt(3)*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/2 + 3*(a + b*x**2)**(sympy.S(1)/3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_658():
    f = (a + b*x**2)**(sympy.S(1)/3)/x**3
    F = -(a + b*x**2)**(sympy.S(1)/3)/(2*x**2) - b*log(x)/(6*a**(sympy.S(2)/3)) + b*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(4*a**(sympy.S(2)/3)) - sqrt(3)*b*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(6*a**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_659():
    f = (a + b*x**2)**(sympy.S(1)/3)/x**5
    F = -(a + b*x**2)**(sympy.S(1)/3)/(4*x**4) - b*(a + b*x**2)**(sympy.S(1)/3)/(12*a*x**2) + b**2*log(x)/(18*a**(sympy.S(5)/3)) - b**2*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(12*a**(sympy.S(5)/3)) + sqrt(3)*b**2*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(18*a**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_660():
    f = x**4*(a + b*x**2)**(sympy.S(1)/3)
    F = -54*3**(sympy.S(3)/4)*a**3*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(935*b**3*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 54*a**2*x*(a + b*x**2)**(sympy.S(1)/3)/(935*b**2) + 6*a*x**3*(a + b*x**2)**(sympy.S(1)/3)/(187*b) + 3*x**5*(a + b*x**2)**(sympy.S(1)/3)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_661():
    f = x**2*(a + b*x**2)**(sympy.S(1)/3)
    F = 6*3**(sympy.S(3)/4)*a**2*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(55*b**2*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + 6*a*x*(a + b*x**2)**(sympy.S(1)/3)/(55*b) + 3*x**3*(a + b*x**2)**(sympy.S(1)/3)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_662():
    f = (a + b*x**2)**(sympy.S(1)/3)
    F = -2*3**(sympy.S(3)/4)*a*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(5*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + 3*x*(a + b*x**2)**(sympy.S(1)/3)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_663():
    f = (a + b*x**2)**(sympy.S(1)/3)/x**2
    F = -2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(3*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - (a + b*x**2)**(sympy.S(1)/3)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_664():
    f = (a + b*x**2)**(sympy.S(1)/3)/x**4
    F = -(a + b*x**2)**(sympy.S(1)/3)/(3*x**3) + 2*3**(sympy.S(3)/4)*b*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(27*a*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 2*b*(a + b*x**2)**(sympy.S(1)/3)/(9*a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_665():
    f = x**7*(a + b*x**2)**(sympy.S(2)/3)
    F = -3*a**3*(a + b*x**2)**(sympy.S(5)/3)/(10*b**4) + 9*a**2*(a + b*x**2)**(sympy.S(8)/3)/(16*b**4) - 9*a*(a + b*x**2)**(sympy.S(11)/3)/(22*b**4) + 3*(a + b*x**2)**(sympy.S(14)/3)/(28*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_666():
    f = x**5*(a + b*x**2)**(sympy.S(2)/3)
    F = 3*a**2*(a + b*x**2)**(sympy.S(5)/3)/(10*b**3) - 3*a*(a + b*x**2)**(sympy.S(8)/3)/(8*b**3) + 3*(a + b*x**2)**(sympy.S(11)/3)/(22*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_667():
    f = x**3*(a + b*x**2)**(sympy.S(2)/3)
    F = -3*a*(a + b*x**2)**(sympy.S(5)/3)/(10*b**2) + 3*(a + b*x**2)**(sympy.S(8)/3)/(16*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_668():
    f = x*(a + b*x**2)**(sympy.S(2)/3)
    F = 3*(a + b*x**2)**(sympy.S(5)/3)/(10*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_669():
    f = (a + b*x**2)**(sympy.S(2)/3)/x
    F = -a**(sympy.S(2)/3)*log(x)/2 + 3*a**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/4 + sqrt(3)*a**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/2 + 3*(a + b*x**2)**(sympy.S(2)/3)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_670():
    f = (a + b*x**2)**(sympy.S(2)/3)/x**3
    F = -(a + b*x**2)**(sympy.S(2)/3)/(2*x**2) - b*log(x)/(3*a**(sympy.S(1)/3)) + b*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(2*a**(sympy.S(1)/3)) + sqrt(3)*b*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_671():
    f = (a + b*x**2)**(sympy.S(2)/3)/x**5
    F = -(a + b*x**2)**(sympy.S(2)/3)/(4*x**4) - b*(a + b*x**2)**(sympy.S(2)/3)/(6*a*x**2) + b**2*log(x)/(18*a**(sympy.S(4)/3)) - b**2*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(12*a**(sympy.S(4)/3)) - sqrt(3)*b**2*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(18*a**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_672():
    f = x**4*(a + b*x**2)**(sympy.S(2)/3)
    F = 162*3**(sympy.S(1)/4)*a**(sympy.S(10)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(1729*b**3*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 108*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(10)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(1729*b**3*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 324*a**3*x/(1729*b**2*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))) - 108*a**2*x*(a + b*x**2)**(sympy.S(2)/3)/(1729*b**2) + 12*a*x**3*(a + b*x**2)**(sympy.S(2)/3)/(247*b) + 3*x**5*(a + b*x**2)**(sympy.S(2)/3)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_673():
    f = x**2*(a + b*x**2)**(sympy.S(2)/3)
    F = -18*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(91*b**2*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + 12*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(91*b**2*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + 36*a**2*x/(91*b*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))) + 12*a*x*(a + b*x**2)**(sympy.S(2)/3)/(91*b) + 3*x**3*(a + b*x**2)**(sympy.S(2)/3)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_674():
    f = (a + b*x**2)**(sympy.S(2)/3)
    F = 6*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(7*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 4*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(7*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 12*a*x/(7*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 7*(a + b*x**2)**(sympy.S(1)/3)) + 3*x*(a + b*x**2)**(sympy.S(2)/3)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_675():
    f = (a + b*x**2)**(sympy.S(2)/3)/x**2
    F = 2*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 4*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(3*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 4*b*x/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3)) - (a + b*x**2)**(sympy.S(2)/3)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_676():
    f = (a + b*x**2)**(sympy.S(2)/3)/x**4
    F = -(a + b*x**2)**(sympy.S(2)/3)/(3*x**3) - 4*b**2*x/(9*a*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))) - 4*b*(a + b*x**2)**(sympy.S(2)/3)/(9*a*x) + 2*3**(sympy.S(1)/4)*b*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(9*a**(sympy.S(2)/3)*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 4*sqrt(2)*3**(sympy.S(3)/4)*b*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(27*a**(sympy.S(2)/3)*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_677():
    f = x**7*(a + b*x**2)**(sympy.S(4)/3)
    F = -3*a**3*(a + b*x**2)**(sympy.S(7)/3)/(14*b**4) + 9*a**2*(a + b*x**2)**(sympy.S(10)/3)/(20*b**4) - 9*a*(a + b*x**2)**(sympy.S(13)/3)/(26*b**4) + 3*(a + b*x**2)**(sympy.S(16)/3)/(32*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_678():
    f = x**5*(a + b*x**2)**(sympy.S(4)/3)
    F = 3*a**2*(a + b*x**2)**(sympy.S(7)/3)/(14*b**3) - 3*a*(a + b*x**2)**(sympy.S(10)/3)/(10*b**3) + 3*(a + b*x**2)**(sympy.S(13)/3)/(26*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_679():
    f = x**3*(a + b*x**2)**(sympy.S(4)/3)
    F = -3*a*(a + b*x**2)**(sympy.S(7)/3)/(14*b**2) + 3*(a + b*x**2)**(sympy.S(10)/3)/(20*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_680():
    f = x*(a + b*x**2)**(sympy.S(4)/3)
    F = 3*(a + b*x**2)**(sympy.S(7)/3)/(14*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_681():
    f = (a + b*x**2)**(sympy.S(4)/3)/x
    F = -a**(sympy.S(4)/3)*log(x)/2 + 3*a**(sympy.S(4)/3)*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/4 - sqrt(3)*a**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/2 + 3*a*(a + b*x**2)**(sympy.S(1)/3)/2 + 3*(a + b*x**2)**(sympy.S(4)/3)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_682():
    f = (a + b*x**2)**(sympy.S(4)/3)/x**3
    F = -2*a**(sympy.S(1)/3)*b*log(x)/3 + a**(sympy.S(1)/3)*b*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3)) - 2*sqrt(3)*a**(sympy.S(1)/3)*b*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/3 + 2*b*(a + b*x**2)**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(4)/3)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_683():
    f = (a + b*x**2)**(sympy.S(4)/3)/x**5
    F = -b*(a + b*x**2)**(sympy.S(1)/3)/(3*x**2) - (a + b*x**2)**(sympy.S(4)/3)/(4*x**4) - b**2*log(x)/(9*a**(sympy.S(2)/3)) + b**2*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(6*a**(sympy.S(2)/3)) - sqrt(3)*b**2*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_684():
    f = x**4*(a + b*x**2)**(sympy.S(4)/3)
    F = -432*3**(sympy.S(3)/4)*a**4*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(21505*b**3*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 432*a**3*x*(a + b*x**2)**(sympy.S(1)/3)/(21505*b**2) + 48*a**2*x**3*(a + b*x**2)**(sympy.S(1)/3)/(4301*b) + 24*a*x**5*(a + b*x**2)**(sympy.S(1)/3)/391 + 3*x**5*(a + b*x**2)**(sympy.S(4)/3)/23
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_685():
    f = x**2*(a + b*x**2)**(sympy.S(4)/3)
    F = 48*3**(sympy.S(3)/4)*a**3*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(935*b**2*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + 48*a**2*x*(a + b*x**2)**(sympy.S(1)/3)/(935*b) + 24*a*x**3*(a + b*x**2)**(sympy.S(1)/3)/187 + 3*x**3*(a + b*x**2)**(sympy.S(4)/3)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_686():
    f = (a + b*x**2)**(sympy.S(4)/3)
    F = -16*3**(sympy.S(3)/4)*a**2*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(55*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + 24*a*x*(a + b*x**2)**(sympy.S(1)/3)/55 + 3*x*(a + b*x**2)**(sympy.S(4)/3)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_687():
    f = (a + b*x**2)**(sympy.S(4)/3)/x**2
    F = -16*3**(sympy.S(3)/4)*a*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(15*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + 8*b*x*(a + b*x**2)**(sympy.S(1)/3)/5 - (a + b*x**2)**(sympy.S(4)/3)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_688():
    f = (a + b*x**2)**(sympy.S(4)/3)/x**4
    F = -16*3**(sympy.S(3)/4)*b*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(27*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 8*b*(a + b*x**2)**(sympy.S(1)/3)/(9*x) - (a + b*x**2)**(sympy.S(4)/3)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_689():
    f = x*(x**2 - 1)**(sympy.S(7)/3)
    F = 3*(x**2 - 1)**(sympy.S(10)/3)/20
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_690():
    f = x**7/(a + b*x**2)**(sympy.S(1)/3)
    F = -3*a**3*(a + b*x**2)**(sympy.S(2)/3)/(4*b**4) + 9*a**2*(a + b*x**2)**(sympy.S(5)/3)/(10*b**4) - 9*a*(a + b*x**2)**(sympy.S(8)/3)/(16*b**4) + 3*(a + b*x**2)**(sympy.S(11)/3)/(22*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_691():
    f = x**5/(a + b*x**2)**(sympy.S(1)/3)
    F = 3*a**2*(a + b*x**2)**(sympy.S(2)/3)/(4*b**3) - 3*a*(a + b*x**2)**(sympy.S(5)/3)/(5*b**3) + 3*(a + b*x**2)**(sympy.S(8)/3)/(16*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_692():
    f = x**3/(a + b*x**2)**(sympy.S(1)/3)
    F = -3*a*(a + b*x**2)**(sympy.S(2)/3)/(4*b**2) + 3*(a + b*x**2)**(sympy.S(5)/3)/(10*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_693():
    f = x/(a + b*x**2)**(sympy.S(1)/3)
    F = 3*(a + b*x**2)**(sympy.S(2)/3)/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_694():
    f = 1/(x*(a + b*x**2)**(sympy.S(1)/3))
    F = -log(x)/(2*a**(sympy.S(1)/3)) + 3*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(4*a**(sympy.S(1)/3)) + sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*a**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_695():
    f = 1/(x**3*(a + b*x**2)**(sympy.S(1)/3))
    F = -(a + b*x**2)**(sympy.S(2)/3)/(2*a*x**2) + b*log(x)/(6*a**(sympy.S(4)/3)) - b*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(4*a**(sympy.S(4)/3)) - sqrt(3)*b*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(6*a**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_696():
    f = 1/(x**5*(a + b*x**2)**(sympy.S(1)/3))
    F = -(a + b*x**2)**(sympy.S(2)/3)/(4*a*x**4) + b*(a + b*x**2)**(sympy.S(2)/3)/(3*a**2*x**2) - b**2*log(x)/(9*a**(sympy.S(7)/3)) + b**2*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(6*a**(sympy.S(7)/3)) + sqrt(3)*b**2*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_697():
    f = x**4/(a + b*x**2)**(sympy.S(1)/3)
    F = 81*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(182*b**3*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 27*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(91*b**3*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 81*a**2*x/(91*b**2*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))) - 27*a*x*(a + b*x**2)**(sympy.S(2)/3)/(91*b**2) + 3*x**3*(a + b*x**2)**(sympy.S(2)/3)/(13*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_698():
    f = x**2/(a + b*x**2)**(sympy.S(1)/3)
    F = -9*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(14*b**2*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + 3*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(7*b**2*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + 9*a*x/(7*b*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))) + 3*x*(a + b*x**2)**(sympy.S(2)/3)/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_699():
    f = (a + b*x**2)**(sympy.S(-1)/3)
    F = 3*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(2*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 3*x/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_700():
    f = 1/(x**2*(a + b*x**2)**(sympy.S(1)/3))
    F = -b*x/(a*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))) - (a + b*x**2)**(sympy.S(2)/3)/(a*x) + 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(2*a**(sympy.S(2)/3)*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(3*a**(sympy.S(2)/3)*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_701():
    f = 1/(x**4*(a + b*x**2)**(sympy.S(1)/3))
    F = -(a + b*x**2)**(sympy.S(2)/3)/(3*a*x**3) + 5*b**2*x/(9*a**2*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))) + 5*b*(a + b*x**2)**(sympy.S(2)/3)/(9*a**2*x) - 5*3**(sympy.S(1)/4)*b*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(18*a**(sympy.S(5)/3)*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + 5*sqrt(2)*3**(sympy.S(3)/4)*b*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(27*a**(sympy.S(5)/3)*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_702():
    f = x**7/(a + b*x**2)**(sympy.S(2)/3)
    F = -3*a**3*(a + b*x**2)**(sympy.S(1)/3)/(2*b**4) + 9*a**2*(a + b*x**2)**(sympy.S(4)/3)/(8*b**4) - 9*a*(a + b*x**2)**(sympy.S(7)/3)/(14*b**4) + 3*(a + b*x**2)**(sympy.S(10)/3)/(20*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_703():
    f = x**5/(a + b*x**2)**(sympy.S(2)/3)
    F = 3*a**2*(a + b*x**2)**(sympy.S(1)/3)/(2*b**3) - 3*a*(a + b*x**2)**(sympy.S(4)/3)/(4*b**3) + 3*(a + b*x**2)**(sympy.S(7)/3)/(14*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_704():
    f = x**3/(a + b*x**2)**(sympy.S(2)/3)
    F = -3*a*(a + b*x**2)**(sympy.S(1)/3)/(2*b**2) + 3*(a + b*x**2)**(sympy.S(4)/3)/(8*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_705():
    f = x/(a + b*x**2)**(sympy.S(2)/3)
    F = 3*(a + b*x**2)**(sympy.S(1)/3)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_706():
    f = 1/(x*(a + b*x**2)**(sympy.S(2)/3))
    F = -log(x)/(2*a**(sympy.S(2)/3)) + 3*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(4*a**(sympy.S(2)/3)) - sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*a**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_707():
    f = 1/(x**3*(a + b*x**2)**(sympy.S(2)/3))
    F = -(a + b*x**2)**(sympy.S(1)/3)/(2*a*x**2) + b*log(x)/(3*a**(sympy.S(5)/3)) - b*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(2*a**(sympy.S(5)/3)) + sqrt(3)*b*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_708():
    f = 1/(x**5*(a + b*x**2)**(sympy.S(2)/3))
    F = -(a + b*x**2)**(sympy.S(1)/3)/(4*a*x**4) + 5*b*(a + b*x**2)**(sympy.S(1)/3)/(12*a**2*x**2) - 5*b**2*log(x)/(18*a**(sympy.S(8)/3)) + 5*b**2*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(12*a**(sympy.S(8)/3)) - 5*sqrt(3)*b**2*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(18*a**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_709():
    f = x**4/(a + b*x**2)**(sympy.S(2)/3)
    F = -27*3**(sympy.S(3)/4)*a**2*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(55*b**3*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 27*a*x*(a + b*x**2)**(sympy.S(1)/3)/(55*b**2) + 3*x**3*(a + b*x**2)**(sympy.S(1)/3)/(11*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_710():
    f = x**2/(a + b*x**2)**(sympy.S(2)/3)
    F = 3*3**(sympy.S(3)/4)*a*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(5*b**2*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + 3*x*(a + b*x**2)**(sympy.S(1)/3)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_711():
    f = (a + b*x**2)**(sympy.S(-2)/3)
    F = -3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_712():
    f = 1/(x**2*(a + b*x**2)**(sympy.S(2)/3))
    F = 3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(3*a*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - (a + b*x**2)**(sympy.S(1)/3)/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_713():
    f = 1/(x**4*(a + b*x**2)**(sympy.S(2)/3))
    F = -(a + b*x**2)**(sympy.S(1)/3)/(3*a*x**3) - 7*3**(sympy.S(3)/4)*b*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(27*a**2*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + 7*b*(a + b*x**2)**(sympy.S(1)/3)/(9*a**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_714():
    f = x**7/(a + b*x**2)**(sympy.S(4)/3)
    F = 3*a**3/(2*b**4*(a + b*x**2)**(sympy.S(1)/3)) + 9*a**2*(a + b*x**2)**(sympy.S(2)/3)/(4*b**4) - 9*a*(a + b*x**2)**(sympy.S(5)/3)/(10*b**4) + 3*(a + b*x**2)**(sympy.S(8)/3)/(16*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_715():
    f = x**5/(a + b*x**2)**(sympy.S(4)/3)
    F = -3*a**2/(2*b**3*(a + b*x**2)**(sympy.S(1)/3)) - 3*a*(a + b*x**2)**(sympy.S(2)/3)/(2*b**3) + 3*(a + b*x**2)**(sympy.S(5)/3)/(10*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_716():
    f = x**3/(a + b*x**2)**(sympy.S(4)/3)
    F = 3*a/(2*b**2*(a + b*x**2)**(sympy.S(1)/3)) + 3*(a + b*x**2)**(sympy.S(2)/3)/(4*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_717():
    f = x/(a + b*x**2)**(sympy.S(4)/3)
    F = -3/(2*b*(a + b*x**2)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_718():
    f = 1/(x*(a + b*x**2)**(sympy.S(4)/3))
    F = 3/(2*a*(a + b*x**2)**(sympy.S(1)/3)) - log(x)/(2*a**(sympy.S(4)/3)) + 3*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(4*a**(sympy.S(4)/3)) + sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*a**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_719():
    f = 1/(x**3*(a + b*x**2)**(sympy.S(4)/3))
    F = -1/(2*a*x**2*(a + b*x**2)**(sympy.S(1)/3)) - 2*b/(a**2*(a + b*x**2)**(sympy.S(1)/3)) + 2*b*log(x)/(3*a**(sympy.S(7)/3)) - b*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/a**(sympy.S(7)/3) - 2*sqrt(3)*b*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_720():
    f = 1/(x**5*(a + b*x**2)**(sympy.S(4)/3))
    F = -1/(4*a*x**4*(a + b*x**2)**(sympy.S(1)/3)) + 7*b/(12*a**2*x**2*(a + b*x**2)**(sympy.S(1)/3)) + 7*b**2/(3*a**3*(a + b*x**2)**(sympy.S(1)/3)) - 7*b**2*log(x)/(9*a**(sympy.S(10)/3)) + 7*b**2*log(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(6*a**(sympy.S(10)/3)) + 7*sqrt(3)*b**2*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(a + b*x**2)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_721():
    f = x**4/(a + b*x**2)**(sympy.S(4)/3)
    F = -81*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(28*b**3*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + 27*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(14*b**3*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + 81*a*x/(14*b**2*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))) - 3*x**3/(2*b*(a + b*x**2)**(sympy.S(1)/3)) + 27*x*(a + b*x**2)**(sympy.S(2)/3)/(14*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_722():
    f = x**2/(a + b*x**2)**(sympy.S(4)/3)
    F = 9*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(4*b**2*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 3*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(2*b**2*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - 9*x/(2*b*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))) - 3*x/(2*b*(a + b*x**2)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_723():
    f = (a + b*x**2)**(sympy.S(-4)/3)
    F = 3*x/(2*a*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))) + 3*x/(2*a*(a + b*x**2)**(sympy.S(1)/3)) - 3*3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(4*a**(sympy.S(2)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(2*a**(sympy.S(2)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_724():
    f = 1/(x**2*(a + b*x**2)**(sympy.S(4)/3))
    F = 3/(2*a*x*(a + b*x**2)**(sympy.S(1)/3)) - 5*b*x/(2*a**2*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))) - 5*(a + b*x**2)**(sympy.S(2)/3)/(2*a**2*x) + 5*3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(4*a**(sympy.S(5)/3)*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) - sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*(5*a**(sympy.S(1)/3) - 5*(a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(6*a**(sympy.S(5)/3)*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_725():
    f = 1/(x**4*(a + b*x**2)**(sympy.S(4)/3))
    F = 3/(2*a*x**3*(a + b*x**2)**(sympy.S(1)/3)) - 11*(a + b*x**2)**(sympy.S(2)/3)/(6*a**2*x**3) + 55*b**2*x/(18*a**3*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))) + 55*b*(a + b*x**2)**(sympy.S(2)/3)/(18*a**3*x) - 55*3**(sympy.S(1)/4)*b*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(36*a**(sympy.S(8)/3)*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)) + 55*sqrt(2)*3**(sympy.S(3)/4)*b*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3) + (a + b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(54*a**(sympy.S(8)/3)*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a + b*x**2)**(sympy.S(1)/3))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_726():
    f = (c*x)**(sympy.S(13)/3)*(a + b*x**2)**(sympy.S(1)/3)
    F = -5*a**3*c**(sympy.S(13)/3)*log(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(162*b**(sympy.S(8)/3)) + 5*a**3*c**(sympy.S(13)/3)*log(b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(324*b**(sympy.S(8)/3)) - 5*sqrt(3)*a**3*c**(sympy.S(13)/3)*atan(sqrt(3)*(2*b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(3*c**(sympy.S(2)/3)))/(162*b**(sympy.S(8)/3)) - 5*a**2*c**3*(c*x)**(sympy.S(4)/3)*(a + b*x**2)**(sympy.S(1)/3)/(108*b**2) + a*c*(c*x)**(sympy.S(10)/3)*(a + b*x**2)**(sympy.S(1)/3)/(36*b) + (c*x)**(sympy.S(16)/3)*(a + b*x**2)**(sympy.S(1)/3)/(6*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_727():
    f = (c*x)**(sympy.S(7)/3)*(a + b*x**2)**(sympy.S(1)/3)
    F = a**2*c**(sympy.S(7)/3)*log(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(18*b**(sympy.S(5)/3)) - a**2*c**(sympy.S(7)/3)*log(b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(36*b**(sympy.S(5)/3)) + sqrt(3)*a**2*c**(sympy.S(7)/3)*atan(sqrt(3)*(2*b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(3*c**(sympy.S(2)/3)))/(18*b**(sympy.S(5)/3)) + a*c*(c*x)**(sympy.S(4)/3)*(a + b*x**2)**(sympy.S(1)/3)/(12*b) + (c*x)**(sympy.S(10)/3)*(a + b*x**2)**(sympy.S(1)/3)/(4*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_728():
    f = (c*x)**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3)
    F = -a*c**(sympy.S(1)/3)*log(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(6*b**(sympy.S(2)/3)) + a*c**(sympy.S(1)/3)*log(b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(12*b**(sympy.S(2)/3)) - sqrt(3)*a*c**(sympy.S(1)/3)*atan(sqrt(3)*(2*b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(3*c**(sympy.S(2)/3)))/(6*b**(sympy.S(2)/3)) + (c*x)**(sympy.S(4)/3)*(a + b*x**2)**(sympy.S(1)/3)/(2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_729():
    f = (a + b*x**2)**(sympy.S(1)/3)/(c*x)**(sympy.S(5)/3)
    F = -b**(sympy.S(1)/3)*log(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(2*c**(sympy.S(5)/3)) + b**(sympy.S(1)/3)*log(b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(4*c**(sympy.S(5)/3)) - sqrt(3)*b**(sympy.S(1)/3)*atan(sqrt(3)*(2*b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(3*c**(sympy.S(2)/3)))/(2*c**(sympy.S(5)/3)) - 3*(a + b*x**2)**(sympy.S(1)/3)/(2*c*(c*x)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_730():
    f = (a + b*x**2)**(sympy.S(1)/3)/(c*x)**(sympy.S(11)/3)
    F = -3*(a + b*x**2)**(sympy.S(4)/3)/(8*a*c*(c*x)**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_731():
    f = (a + b*x**2)**(sympy.S(1)/3)/(c*x)**(sympy.S(17)/3)
    F = -3*(a + b*x**2)**(sympy.S(4)/3)/(8*a*c*(c*x)**(sympy.S(14)/3)) + 9*(a + b*x**2)**(sympy.S(7)/3)/(56*a**2*c*(c*x)**(sympy.S(14)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_732():
    f = (a + b*x**2)**(sympy.S(1)/3)/(c*x)**(sympy.S(23)/3)
    F = -3*(a + b*x**2)**(sympy.S(4)/3)/(8*a*c*(c*x)**(sympy.S(20)/3)) + 9*(a + b*x**2)**(sympy.S(7)/3)/(28*a**2*c*(c*x)**(sympy.S(20)/3)) - 27*(a + b*x**2)**(sympy.S(10)/3)/(280*a**3*c*(c*x)**(sympy.S(20)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_733():
    f = (a + b*x**2)**(sympy.S(1)/3)/(c*x)**(sympy.S(29)/3)
    F = -3*(a + b*x**2)**(sympy.S(4)/3)/(8*a*c*(c*x)**(sympy.S(26)/3)) + 27*(a + b*x**2)**(sympy.S(7)/3)/(56*a**2*c*(c*x)**(sympy.S(26)/3)) - 81*(a + b*x**2)**(sympy.S(10)/3)/(280*a**3*c*(c*x)**(sympy.S(26)/3)) + 243*(a + b*x**2)**(sympy.S(13)/3)/(3640*a**4*c*(c*x)**(sympy.S(26)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_734():
    f = (c*x)**(sympy.S(10)/3)*(a + b*x**2)**(sympy.S(1)/3)
    F = 7*3**(sympy.S(3)/4)*a**2*c**(sympy.S(7)/3)*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(405*b**2*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2))) - 14*a**2*c**3*(c*x)**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3)/(135*b**2) + 2*a*c*(c*x)**(sympy.S(7)/3)*(a + b*x**2)**(sympy.S(1)/3)/(45*b) + (c*x)**(sympy.S(13)/3)*(a + b*x**2)**(sympy.S(1)/3)/(5*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_735():
    f = (c*x)**(sympy.S(4)/3)*(a + b*x**2)**(sympy.S(1)/3)
    F = -3**(sympy.S(3)/4)*a*c**(sympy.S(1)/3)*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(27*b*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2))) + 2*a*c*(c*x)**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3)/(9*b) + (c*x)**(sympy.S(7)/3)*(a + b*x**2)**(sympy.S(1)/3)/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_736():
    f = (a + b*x**2)**(sympy.S(1)/3)/(c*x)**(sympy.S(2)/3)
    F = (c*x)**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3)/c + 3**(sympy.S(3)/4)*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(3*c**(sympy.S(5)/3)*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_737():
    f = (a + b*x**2)**(sympy.S(1)/3)/(c*x)**(sympy.S(8)/3)
    F = -3*(a + b*x**2)**(sympy.S(1)/3)/(5*c*(c*x)**(sympy.S(5)/3)) + 3**(sympy.S(3)/4)*b*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(5*a*c**(sympy.S(11)/3)*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_738():
    f = (a + b*x**2)**(sympy.S(1)/3)/(c*x)**(sympy.S(14)/3)
    F = -3*(a + b*x**2)**(sympy.S(1)/3)/(11*c*(c*x)**(sympy.S(11)/3)) - 6*b*(a + b*x**2)**(sympy.S(1)/3)/(55*a*c**3*(c*x)**(sympy.S(5)/3)) - 3*3**(sympy.S(3)/4)*b**2*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(55*a**2*c**(sympy.S(17)/3)*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_739():
    f = (c*x)**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(1)/3)
    F = 3*(c*x)**(sympy.S(5)/3)*(a + b*x**2)**(sympy.S(1)/3)*hyper((sympy.S(-1)/3, sympy.S(5)/6), (sympy.S(11)/6,), -b*x**2/a)/(5*c*(1 + b*x**2/a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_740():
    f = (a + b*x**2)**(sympy.S(1)/3)/(c*x)**(sympy.S(1)/3)
    F = 3*(c*x)**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(1)/3)*hyper((sympy.S(-1)/3, sympy.S(1)/3), (sympy.S(4)/3,), -b*x**2/a)/(2*c*(1 + b*x**2/a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_741():
    f = (a + b*x**2)**(sympy.S(1)/3)/(c*x)**(sympy.S(4)/3)
    F = -3*(a + b*x**2)**(sympy.S(1)/3)*hyper((sympy.S(-1)/3, sympy.S(-1)/6), (sympy.S(5)/6,), -b*x**2/a)/(c*(c*x)**(sympy.S(1)/3)*(1 + b*x**2/a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_742():
    f = (c*x)**(sympy.S(13)/3)*(a + b*x**2)**(sympy.S(4)/3)
    F = -5*a**4*c**(sympy.S(13)/3)*log(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(486*b**(sympy.S(8)/3)) + 5*a**4*c**(sympy.S(13)/3)*log(b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(972*b**(sympy.S(8)/3)) - 5*sqrt(3)*a**4*c**(sympy.S(13)/3)*atan(sqrt(3)*(2*b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(3*c**(sympy.S(2)/3)))/(486*b**(sympy.S(8)/3)) - 5*a**3*c**3*(c*x)**(sympy.S(4)/3)*(a + b*x**2)**(sympy.S(1)/3)/(324*b**2) + a**2*c*(c*x)**(sympy.S(10)/3)*(a + b*x**2)**(sympy.S(1)/3)/(108*b) + a*(c*x)**(sympy.S(16)/3)*(a + b*x**2)**(sympy.S(1)/3)/(18*c) + (c*x)**(sympy.S(16)/3)*(a + b*x**2)**(sympy.S(4)/3)/(8*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_743():
    f = (c*x)**(sympy.S(7)/3)*(a + b*x**2)**(sympy.S(4)/3)
    F = 2*a**3*c**(sympy.S(7)/3)*log(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(81*b**(sympy.S(5)/3)) - a**3*c**(sympy.S(7)/3)*log(b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(81*b**(sympy.S(5)/3)) + 2*sqrt(3)*a**3*c**(sympy.S(7)/3)*atan(sqrt(3)*(2*b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(3*c**(sympy.S(2)/3)))/(81*b**(sympy.S(5)/3)) + a**2*c*(c*x)**(sympy.S(4)/3)*(a + b*x**2)**(sympy.S(1)/3)/(27*b) + a*(c*x)**(sympy.S(10)/3)*(a + b*x**2)**(sympy.S(1)/3)/(9*c) + (c*x)**(sympy.S(10)/3)*(a + b*x**2)**(sympy.S(4)/3)/(6*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_744():
    f = (c*x)**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(4)/3)
    F = -a**2*c**(sympy.S(1)/3)*log(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(9*b**(sympy.S(2)/3)) + a**2*c**(sympy.S(1)/3)*log(b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(18*b**(sympy.S(2)/3)) - sqrt(3)*a**2*c**(sympy.S(1)/3)*atan(sqrt(3)*(2*b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(3*c**(sympy.S(2)/3)))/(9*b**(sympy.S(2)/3)) + a*(c*x)**(sympy.S(4)/3)*(a + b*x**2)**(sympy.S(1)/3)/(3*c) + (c*x)**(sympy.S(4)/3)*(a + b*x**2)**(sympy.S(4)/3)/(4*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_745():
    f = (a + b*x**2)**(sympy.S(4)/3)/(c*x)**(sympy.S(5)/3)
    F = -2*a*b**(sympy.S(1)/3)*log(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(3*c**(sympy.S(5)/3)) + a*b**(sympy.S(1)/3)*log(b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(3*c**(sympy.S(5)/3)) - 2*sqrt(3)*a*b**(sympy.S(1)/3)*atan(sqrt(3)*(2*b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(3*c**(sympy.S(2)/3)))/(3*c**(sympy.S(5)/3)) + 2*b*(c*x)**(sympy.S(4)/3)*(a + b*x**2)**(sympy.S(1)/3)/c**3 - 3*(a + b*x**2)**(sympy.S(4)/3)/(2*c*(c*x)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_746():
    f = (a + b*x**2)**(sympy.S(4)/3)/(c*x)**(sympy.S(11)/3)
    F = -b**(sympy.S(4)/3)*log(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(2*c**(sympy.S(11)/3)) + b**(sympy.S(4)/3)*log(b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(4*c**(sympy.S(11)/3)) - sqrt(3)*b**(sympy.S(4)/3)*atan(sqrt(3)*(2*b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(3*c**(sympy.S(2)/3)))/(2*c**(sympy.S(11)/3)) - 3*b*(a + b*x**2)**(sympy.S(1)/3)/(2*c**3*(c*x)**(sympy.S(2)/3)) - 3*(a + b*x**2)**(sympy.S(4)/3)/(8*c*(c*x)**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_747():
    f = (a + b*x**2)**(sympy.S(4)/3)/(c*x)**(sympy.S(17)/3)
    F = -3*(a + b*x**2)**(sympy.S(7)/3)/(14*a*c*(c*x)**(sympy.S(14)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_748():
    f = (a + b*x**2)**(sympy.S(4)/3)/(c*x)**(sympy.S(23)/3)
    F = -3*(a + b*x**2)**(sympy.S(7)/3)/(14*a*c*(c*x)**(sympy.S(20)/3)) + 9*(a + b*x**2)**(sympy.S(10)/3)/(140*a**2*c*(c*x)**(sympy.S(20)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_749():
    f = (a + b*x**2)**(sympy.S(4)/3)/(c*x)**(sympy.S(29)/3)
    F = -3*(a + b*x**2)**(sympy.S(7)/3)/(14*a*c*(c*x)**(sympy.S(26)/3)) + 9*(a + b*x**2)**(sympy.S(10)/3)/(70*a**2*c*(c*x)**(sympy.S(26)/3)) - 27*(a + b*x**2)**(sympy.S(13)/3)/(910*a**3*c*(c*x)**(sympy.S(26)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_750():
    f = (c*x)**(sympy.S(10)/3)*(a + b*x**2)**(sympy.S(4)/3)
    F = 8*3**(sympy.S(3)/4)*a**3*c**(sympy.S(7)/3)*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(1215*b**2*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2))) - 16*a**3*c**3*(c*x)**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3)/(405*b**2) + 16*a**2*c*(c*x)**(sympy.S(7)/3)*(a + b*x**2)**(sympy.S(1)/3)/(945*b) + 8*a*(c*x)**(sympy.S(13)/3)*(a + b*x**2)**(sympy.S(1)/3)/(105*c) + (c*x)**(sympy.S(13)/3)*(a + b*x**2)**(sympy.S(4)/3)/(7*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_751():
    f = (c*x)**(sympy.S(4)/3)*(a + b*x**2)**(sympy.S(4)/3)
    F = -8*3**(sympy.S(3)/4)*a**2*c**(sympy.S(1)/3)*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(405*b*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2))) + 16*a**2*c*(c*x)**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3)/(135*b) + 8*a*(c*x)**(sympy.S(7)/3)*(a + b*x**2)**(sympy.S(1)/3)/(45*c) + (c*x)**(sympy.S(7)/3)*(a + b*x**2)**(sympy.S(4)/3)/(5*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_752():
    f = (a + b*x**2)**(sympy.S(4)/3)/(c*x)**(sympy.S(2)/3)
    F = 8*a*(c*x)**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3)/(9*c) + 8*3**(sympy.S(3)/4)*a*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(27*c**(sympy.S(5)/3)*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2))) + (c*x)**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(4)/3)/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_753():
    f = (a + b*x**2)**(sympy.S(4)/3)/(c*x)**(sympy.S(8)/3)
    F = 8*b*(c*x)**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3)/(5*c**3) + 8*3**(sympy.S(3)/4)*b*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(15*c**(sympy.S(11)/3)*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2))) - 3*(a + b*x**2)**(sympy.S(4)/3)/(5*c*(c*x)**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_754():
    f = (a + b*x**2)**(sympy.S(4)/3)/(c*x)**(sympy.S(14)/3)
    F = -24*b*(a + b*x**2)**(sympy.S(1)/3)/(55*c**3*(c*x)**(sympy.S(5)/3)) - 3*(a + b*x**2)**(sympy.S(4)/3)/(11*c*(c*x)**(sympy.S(11)/3)) + 8*3**(sympy.S(3)/4)*b**2*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(55*a*c**(sympy.S(17)/3)*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_755():
    f = (a + b*x**2)**(sympy.S(4)/3)/(c*x)**(sympy.S(20)/3)
    F = -24*b*(a + b*x**2)**(sympy.S(1)/3)/(187*c**3*(c*x)**(sympy.S(11)/3)) - 3*(a + b*x**2)**(sympy.S(4)/3)/(17*c*(c*x)**(sympy.S(17)/3)) - 48*b**2*(a + b*x**2)**(sympy.S(1)/3)/(935*a*c**5*(c*x)**(sympy.S(5)/3)) - 24*3**(sympy.S(3)/4)*b**3*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(935*a**2*c**(sympy.S(23)/3)*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_756():
    f = (c*x)**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(4)/3)
    F = 3*a*(c*x)**(sympy.S(5)/3)*(a + b*x**2)**(sympy.S(1)/3)*hyper((sympy.S(-4)/3, sympy.S(5)/6), (sympy.S(11)/6,), -b*x**2/a)/(5*c*(1 + b*x**2/a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_757():
    f = (a + b*x**2)**(sympy.S(4)/3)/(c*x)**(sympy.S(1)/3)
    F = 3*a*(c*x)**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(1)/3)*hyper((sympy.S(-4)/3, sympy.S(1)/3), (sympy.S(4)/3,), -b*x**2/a)/(2*c*(1 + b*x**2/a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_758():
    f = (a + b*x**2)**(sympy.S(4)/3)/(c*x)**(sympy.S(4)/3)
    F = -3*a*(a + b*x**2)**(sympy.S(1)/3)*hyper((sympy.S(-4)/3, sympy.S(-1)/6), (sympy.S(5)/6,), -b*x**2/a)/(c*(c*x)**(sympy.S(1)/3)*(1 + b*x**2/a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_759():
    f = (c*x)**(sympy.S(19)/3)/(a + b*x**2)**(sympy.S(2)/3)
    F = 20*a**3*c**(sympy.S(19)/3)*log(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(81*b**(sympy.S(11)/3)) - 10*a**3*c**(sympy.S(19)/3)*log(b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(81*b**(sympy.S(11)/3)) + 20*sqrt(3)*a**3*c**(sympy.S(19)/3)*atan(sqrt(3)*(2*b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(3*c**(sympy.S(2)/3)))/(81*b**(sympy.S(11)/3)) + 10*a**2*c**5*(c*x)**(sympy.S(4)/3)*(a + b*x**2)**(sympy.S(1)/3)/(27*b**3) - 2*a*c**3*(c*x)**(sympy.S(10)/3)*(a + b*x**2)**(sympy.S(1)/3)/(9*b**2) + c*(c*x)**(sympy.S(16)/3)*(a + b*x**2)**(sympy.S(1)/3)/(6*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_760():
    f = (c*x)**(sympy.S(13)/3)/(a + b*x**2)**(sympy.S(2)/3)
    F = -5*a**2*c**(sympy.S(13)/3)*log(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(18*b**(sympy.S(8)/3)) + 5*a**2*c**(sympy.S(13)/3)*log(b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(36*b**(sympy.S(8)/3)) - 5*sqrt(3)*a**2*c**(sympy.S(13)/3)*atan(sqrt(3)*(2*b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(3*c**(sympy.S(2)/3)))/(18*b**(sympy.S(8)/3)) - 5*a*c**3*(c*x)**(sympy.S(4)/3)*(a + b*x**2)**(sympy.S(1)/3)/(12*b**2) + c*(c*x)**(sympy.S(10)/3)*(a + b*x**2)**(sympy.S(1)/3)/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_761():
    f = (c*x)**(sympy.S(7)/3)/(a + b*x**2)**(sympy.S(2)/3)
    F = a*c**(sympy.S(7)/3)*log(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(3*b**(sympy.S(5)/3)) - a*c**(sympy.S(7)/3)*log(b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(6*b**(sympy.S(5)/3)) + sqrt(3)*a*c**(sympy.S(7)/3)*atan(sqrt(3)*(2*b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(3*c**(sympy.S(2)/3)))/(3*b**(sympy.S(5)/3)) + c*(c*x)**(sympy.S(4)/3)*(a + b*x**2)**(sympy.S(1)/3)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_762():
    f = (c*x)**(sympy.S(1)/3)/(a + b*x**2)**(sympy.S(2)/3)
    F = -c**(sympy.S(1)/3)*log(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(2*b**(sympy.S(2)/3)) + c**(sympy.S(1)/3)*log(b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(4*b**(sympy.S(2)/3)) - sqrt(3)*c**(sympy.S(1)/3)*atan(sqrt(3)*(2*b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(3*c**(sympy.S(2)/3)))/(2*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_763():
    f = 1/((c*x)**(sympy.S(5)/3)*(a + b*x**2)**(sympy.S(2)/3))
    F = -3*(a + b*x**2)**(sympy.S(1)/3)/(2*a*c*(c*x)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_764():
    f = 1/((c*x)**(sympy.S(11)/3)*(a + b*x**2)**(sympy.S(2)/3))
    F = -3*(a + b*x**2)**(sympy.S(1)/3)/(2*a*c*(c*x)**(sympy.S(8)/3)) + 9*(a + b*x**2)**(sympy.S(4)/3)/(8*a**2*c*(c*x)**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_765():
    f = 1/((c*x)**(sympy.S(17)/3)*(a + b*x**2)**(sympy.S(2)/3))
    F = -3*(a + b*x**2)**(sympy.S(1)/3)/(2*a*c*(c*x)**(sympy.S(14)/3)) + 9*(a + b*x**2)**(sympy.S(4)/3)/(4*a**2*c*(c*x)**(sympy.S(14)/3)) - 27*(a + b*x**2)**(sympy.S(7)/3)/(28*a**3*c*(c*x)**(sympy.S(14)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_766():
    f = 1/((c*x)**(sympy.S(23)/3)*(a + b*x**2)**(sympy.S(2)/3))
    F = -3*(a + b*x**2)**(sympy.S(1)/3)/(2*a*c*(c*x)**(sympy.S(20)/3)) + 27*(a + b*x**2)**(sympy.S(4)/3)/(8*a**2*c*(c*x)**(sympy.S(20)/3)) - 81*(a + b*x**2)**(sympy.S(7)/3)/(28*a**3*c*(c*x)**(sympy.S(20)/3)) + 243*(a + b*x**2)**(sympy.S(10)/3)/(280*a**4*c*(c*x)**(sympy.S(20)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_767():
    f = (c*x)**(sympy.S(10)/3)/(a + b*x**2)**(sympy.S(2)/3)
    F = 7*3**(sympy.S(3)/4)*a*c**(sympy.S(7)/3)*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(54*b**2*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2))) - 7*a*c**3*(c*x)**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3)/(9*b**2) + c*(c*x)**(sympy.S(7)/3)*(a + b*x**2)**(sympy.S(1)/3)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_768():
    f = (c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3)
    F = -3**(sympy.S(3)/4)*c**(sympy.S(1)/3)*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(6*b*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2))) + c*(c*x)**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_769():
    f = 1/((c*x)**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(2)/3))
    F = 3**(sympy.S(3)/4)*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(2*a*c**(sympy.S(5)/3)*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_770():
    f = 1/((c*x)**(sympy.S(8)/3)*(a + b*x**2)**(sympy.S(2)/3))
    F = -3*(a + b*x**2)**(sympy.S(1)/3)/(5*a*c*(c*x)**(sympy.S(5)/3)) - 3*3**(sympy.S(3)/4)*b*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(10*a**2*c**(sympy.S(11)/3)*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_771():
    f = 1/((c*x)**(sympy.S(14)/3)*(a + b*x**2)**(sympy.S(2)/3))
    F = -3*(a + b*x**2)**(sympy.S(1)/3)/(11*a*c*(c*x)**(sympy.S(11)/3)) + 27*b*(a + b*x**2)**(sympy.S(1)/3)/(55*a**2*c**3*(c*x)**(sympy.S(5)/3)) + 27*3**(sympy.S(3)/4)*b**2*(c*x)**(sympy.S(1)/3)*sqrt((b**(sympy.S(2)/3)*(c*x)**(sympy.S(4)/3)/(a + b*x**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(4)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)*(a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))*elliptic_f(acos((-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 - sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))), sqrt(3)/4 + sympy.S.Half)/(110*a**3*c**(sympy.S(17)/3)*sqrt(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))/((a + b*x**2)**(sympy.S(1)/3)*(-b**(sympy.S(1)/3)*(c*x)**(sympy.S(2)/3)*(1 + sqrt(3))/(a + b*x**2)**(sympy.S(1)/3) + c**(sympy.S(2)/3))**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_772():
    f = (c*x)**(sympy.S(2)/3)/(a + b*x**2)**(sympy.S(2)/3)
    F = 3*(c*x)**(sympy.S(5)/3)*(1 + b*x**2/a)**(sympy.S(2)/3)*hyper((sympy.S(2)/3, sympy.S(5)/6), (sympy.S(11)/6,), -b*x**2/a)/(5*c*(a + b*x**2)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_773():
    f = 1/((c*x)**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(2)/3))
    F = 3*(c*x)**(sympy.S(2)/3)*(1 + b*x**2/a)**(sympy.S(2)/3)*hyper((sympy.S(1)/3, sympy.S(2)/3), (sympy.S(4)/3,), -b*x**2/a)/(2*c*(a + b*x**2)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_774():
    f = 1/((c*x)**(sympy.S(4)/3)*(a + b*x**2)**(sympy.S(2)/3))
    F = -3*(1 + b*x**2/a)**(sympy.S(2)/3)*hyper((sympy.S(-1)/6, sympy.S(2)/3), (sympy.S(5)/6,), -b*x**2/a)/(c*(c*x)**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_775():
    f = x**4*(a + b*x**2)**(sympy.S(1)/4)
    F = 8*a**(sympy.S(7)/2)*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(77*b**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/4)) - 4*a**2*x*(a + b*x**2)**(sympy.S(1)/4)/(77*b**2) + 2*a*x**3*(a + b*x**2)**(sympy.S(1)/4)/(77*b) + 2*x**5*(a + b*x**2)**(sympy.S(1)/4)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_776():
    f = x**2*(a + b*x**2)**(sympy.S(1)/4)
    F = -4*a**(sympy.S(5)/2)*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(21*b**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/4)) + 2*a*x*(a + b*x**2)**(sympy.S(1)/4)/(21*b) + 2*x**3*(a + b*x**2)**(sympy.S(1)/4)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_777():
    f = (a + b*x**2)**(sympy.S(1)/4)
    F = 2*a**(sympy.S(3)/2)*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(3*sqrt(b)*(a + b*x**2)**(sympy.S(3)/4)) + 2*x*(a + b*x**2)**(sympy.S(1)/4)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_778():
    f = (a + b*x**2)**(sympy.S(1)/4)/x**2
    F = sqrt(a)*sqrt(b)*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(a + b*x**2)**(sympy.S(3)/4) - (a + b*x**2)**(sympy.S(1)/4)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_779():
    f = (a + b*x**2)**(sympy.S(1)/4)/x**4
    F = -(a + b*x**2)**(sympy.S(1)/4)/(3*x**3) - b*(a + b*x**2)**(sympy.S(1)/4)/(6*a*x) - b**(sympy.S(3)/2)*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(6*sqrt(a)*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_780():
    f = (a + b*x**2)**(sympy.S(1)/4)/x**6
    F = -(a + b*x**2)**(sympy.S(1)/4)/(5*x**5) - b*(a + b*x**2)**(sympy.S(1)/4)/(30*a*x**3) + b**2*(a + b*x**2)**(sympy.S(1)/4)/(12*a**2*x) + b**(sympy.S(5)/2)*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(12*a**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_781():
    f = x**4*(a - b*x**2)**(sympy.S(1)/4)
    F = 8*a**(sympy.S(7)/2)*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(77*b**(sympy.S(5)/2)*(a - b*x**2)**(sympy.S(3)/4)) - 4*a**2*x*(a - b*x**2)**(sympy.S(1)/4)/(77*b**2) - 2*a*x**3*(a - b*x**2)**(sympy.S(1)/4)/(77*b) + 2*x**5*(a - b*x**2)**(sympy.S(1)/4)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_782():
    f = x**2*(a - b*x**2)**(sympy.S(1)/4)
    F = 4*a**(sympy.S(5)/2)*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(21*b**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(3)/4)) - 2*a*x*(a - b*x**2)**(sympy.S(1)/4)/(21*b) + 2*x**3*(a - b*x**2)**(sympy.S(1)/4)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_783():
    f = (a - b*x**2)**(sympy.S(1)/4)
    F = 2*a**(sympy.S(3)/2)*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(3*sqrt(b)*(a - b*x**2)**(sympy.S(3)/4)) + 2*x*(a - b*x**2)**(sympy.S(1)/4)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_784():
    f = (a - b*x**2)**(sympy.S(1)/4)/x**2
    F = -sqrt(a)*sqrt(b)*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(a - b*x**2)**(sympy.S(3)/4) - (a - b*x**2)**(sympy.S(1)/4)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_785():
    f = (a - b*x**2)**(sympy.S(1)/4)/x**4
    F = -(a - b*x**2)**(sympy.S(1)/4)/(3*x**3) + b*(a - b*x**2)**(sympy.S(1)/4)/(6*a*x) - b**(sympy.S(3)/2)*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(6*sqrt(a)*(a - b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_786():
    f = (a - b*x**2)**(sympy.S(1)/4)/x**6
    F = -(a - b*x**2)**(sympy.S(1)/4)/(5*x**5) + b*(a - b*x**2)**(sympy.S(1)/4)/(30*a*x**3) + b**2*(a - b*x**2)**(sympy.S(1)/4)/(12*a**2*x) - b**(sympy.S(5)/2)*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(12*a**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_787():
    f = x**4*(a + b*x**2)**(sympy.S(3)/4)
    F = -8*a**(sympy.S(7)/2)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(65*b**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)) + 8*a**3*x/(65*b**2*(a + b*x**2)**(sympy.S(1)/4)) - 4*a**2*x*(a + b*x**2)**(sympy.S(3)/4)/(65*b**2) + 2*a*x**3*(a + b*x**2)**(sympy.S(3)/4)/(39*b) + 2*x**5*(a + b*x**2)**(sympy.S(3)/4)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_788():
    f = x**2*(a + b*x**2)**(sympy.S(3)/4)
    F = 4*a**(sympy.S(5)/2)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(15*b**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)) - 4*a**2*x/(15*b*(a + b*x**2)**(sympy.S(1)/4)) + 2*a*x*(a + b*x**2)**(sympy.S(3)/4)/(15*b) + 2*x**3*(a + b*x**2)**(sympy.S(3)/4)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_789():
    f = (a + b*x**2)**(sympy.S(3)/4)
    F = -6*a**(sympy.S(3)/2)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(5*sqrt(b)*(a + b*x**2)**(sympy.S(1)/4)) + 6*a*x/(5*(a + b*x**2)**(sympy.S(1)/4)) + 2*x*(a + b*x**2)**(sympy.S(3)/4)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_790():
    f = (a + b*x**2)**(sympy.S(3)/4)/x**2
    F = -3*sqrt(a)*sqrt(b)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(a + b*x**2)**(sympy.S(1)/4) + 3*b*x/(a + b*x**2)**(sympy.S(1)/4) - (a + b*x**2)**(sympy.S(3)/4)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_791():
    f = (a + b*x**2)**(sympy.S(3)/4)/x**4
    F = -(a + b*x**2)**(sympy.S(3)/4)/(3*x**3) + b**2*x/(2*a*(a + b*x**2)**(sympy.S(1)/4)) - b*(a + b*x**2)**(sympy.S(3)/4)/(2*a*x) - b**(sympy.S(3)/2)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(2*sqrt(a)*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_792():
    f = (a + b*x**2)**(sympy.S(3)/4)/x**6
    F = -(a + b*x**2)**(sympy.S(3)/4)/(5*x**5) - b*(a + b*x**2)**(sympy.S(3)/4)/(10*a*x**3) - 3*b**3*x/(20*a**2*(a + b*x**2)**(sympy.S(1)/4)) + 3*b**2*(a + b*x**2)**(sympy.S(3)/4)/(20*a**2*x) + 3*b**(sympy.S(5)/2)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(20*a**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_793():
    f = x**4*(a - b*x**2)**(sympy.S(3)/4)
    F = 8*a**(sympy.S(7)/2)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(65*b**(sympy.S(5)/2)*(a - b*x**2)**(sympy.S(1)/4)) - 4*a**2*x*(a - b*x**2)**(sympy.S(3)/4)/(65*b**2) - 2*a*x**3*(a - b*x**2)**(sympy.S(3)/4)/(39*b) + 2*x**5*(a - b*x**2)**(sympy.S(3)/4)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_794():
    f = x**2*(a - b*x**2)**(sympy.S(3)/4)
    F = 4*a**(sympy.S(5)/2)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(15*b**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(1)/4)) - 2*a*x*(a - b*x**2)**(sympy.S(3)/4)/(15*b) + 2*x**3*(a - b*x**2)**(sympy.S(3)/4)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_795():
    f = (a - b*x**2)**(sympy.S(3)/4)
    F = 6*a**(sympy.S(3)/2)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(5*sqrt(b)*(a - b*x**2)**(sympy.S(1)/4)) + 2*x*(a - b*x**2)**(sympy.S(3)/4)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_796():
    f = (a - b*x**2)**(sympy.S(3)/4)/x**2
    F = -3*sqrt(a)*sqrt(b)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(a - b*x**2)**(sympy.S(1)/4) - (a - b*x**2)**(sympy.S(3)/4)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_797():
    f = (a - b*x**2)**(sympy.S(3)/4)/x**4
    F = -(a - b*x**2)**(sympy.S(3)/4)/(3*x**3) + b*(a - b*x**2)**(sympy.S(3)/4)/(2*a*x) + b**(sympy.S(3)/2)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(2*sqrt(a)*(a - b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_798():
    f = (a - b*x**2)**(sympy.S(3)/4)/x**6
    F = -(a - b*x**2)**(sympy.S(3)/4)/(5*x**5) + b*(a - b*x**2)**(sympy.S(3)/4)/(10*a*x**3) + 3*b**2*(a - b*x**2)**(sympy.S(3)/4)/(20*a**2*x) + 3*b**(sympy.S(5)/2)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(20*a**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_799():
    f = (a + b*x**2)**(sympy.S(5)/4)
    F = 10*a**(sympy.S(5)/2)*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(21*sqrt(b)*(a + b*x**2)**(sympy.S(3)/4)) + 10*a*x*(a + b*x**2)**(sympy.S(1)/4)/21 + 2*x*(a + b*x**2)**(sympy.S(5)/4)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_800():
    f = (a - b*x**2)**(sympy.S(5)/4)
    F = 10*a**(sympy.S(5)/2)*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(21*sqrt(b)*(a - b*x**2)**(sympy.S(3)/4)) + 10*a*x*(a - b*x**2)**(sympy.S(1)/4)/21 + 2*x*(a - b*x**2)**(sympy.S(5)/4)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_801():
    f = (a + b*x**2)**(sympy.S(7)/4)
    F = -14*a**(sympy.S(5)/2)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(15*sqrt(b)*(a + b*x**2)**(sympy.S(1)/4)) + 14*a**2*x/(15*(a + b*x**2)**(sympy.S(1)/4)) + 14*a*x*(a + b*x**2)**(sympy.S(3)/4)/45 + 2*x*(a + b*x**2)**(sympy.S(7)/4)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_802():
    f = (a - b*x**2)**(sympy.S(7)/4)
    F = 14*a**(sympy.S(5)/2)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(15*sqrt(b)*(a - b*x**2)**(sympy.S(1)/4)) + 14*a*x*(a - b*x**2)**(sympy.S(3)/4)/45 + 2*x*(a - b*x**2)**(sympy.S(7)/4)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_803():
    f = x**6/(a + b*x**2)**(sympy.S(1)/4)
    F = 16*a**(sympy.S(7)/2)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(39*b**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(1)/4)) - 16*a**3*x/(39*b**3*(a + b*x**2)**(sympy.S(1)/4)) + 8*a**2*x*(a + b*x**2)**(sympy.S(3)/4)/(39*b**3) - 20*a*x**3*(a + b*x**2)**(sympy.S(3)/4)/(117*b**2) + 2*x**5*(a + b*x**2)**(sympy.S(3)/4)/(13*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_804():
    f = x**4/(a + b*x**2)**(sympy.S(1)/4)
    F = -8*a**(sympy.S(5)/2)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(15*b**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)) + 8*a**2*x/(15*b**2*(a + b*x**2)**(sympy.S(1)/4)) - 4*a*x*(a + b*x**2)**(sympy.S(3)/4)/(15*b**2) + 2*x**3*(a + b*x**2)**(sympy.S(3)/4)/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_805():
    f = x**2/(a + b*x**2)**(sympy.S(1)/4)
    F = 4*a**(sympy.S(3)/2)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(5*b**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)) - 4*a*x/(5*b*(a + b*x**2)**(sympy.S(1)/4)) + 2*x*(a + b*x**2)**(sympy.S(3)/4)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_806():
    f = (a + b*x**2)**(sympy.S(-1)/4)
    F = -2*sqrt(a)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(b)*(a + b*x**2)**(sympy.S(1)/4)) + 2*x/(a + b*x**2)**(sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_807():
    f = 1/(x**2*(a + b*x**2)**(sympy.S(1)/4))
    F = b*x/(a*(a + b*x**2)**(sympy.S(1)/4)) - (a + b*x**2)**(sympy.S(3)/4)/(a*x) - sqrt(b)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(a)*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_808():
    f = 1/(x**4*(a + b*x**2)**(sympy.S(1)/4))
    F = -(a + b*x**2)**(sympy.S(3)/4)/(3*a*x**3) - b**2*x/(2*a**2*(a + b*x**2)**(sympy.S(1)/4)) + b*(a + b*x**2)**(sympy.S(3)/4)/(2*a**2*x) + b**(sympy.S(3)/2)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(2*a**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_809():
    f = 1/(x**6*(a + b*x**2)**(sympy.S(1)/4))
    F = -(a + b*x**2)**(sympy.S(3)/4)/(5*a*x**5) + 7*b*(a + b*x**2)**(sympy.S(3)/4)/(30*a**2*x**3) + 7*b**3*x/(20*a**3*(a + b*x**2)**(sympy.S(1)/4)) - 7*b**2*(a + b*x**2)**(sympy.S(3)/4)/(20*a**3*x) - 7*b**(sympy.S(5)/2)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(20*a**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_810():
    f = x**6/(a - b*x**2)**(sympy.S(1)/4)
    F = 16*a**(sympy.S(7)/2)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(39*b**(sympy.S(7)/2)*(a - b*x**2)**(sympy.S(1)/4)) - 8*a**2*x*(a - b*x**2)**(sympy.S(3)/4)/(39*b**3) - 20*a*x**3*(a - b*x**2)**(sympy.S(3)/4)/(117*b**2) - 2*x**5*(a - b*x**2)**(sympy.S(3)/4)/(13*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_811():
    f = x**4/(a - b*x**2)**(sympy.S(1)/4)
    F = 8*a**(sympy.S(5)/2)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(15*b**(sympy.S(5)/2)*(a - b*x**2)**(sympy.S(1)/4)) - 4*a*x*(a - b*x**2)**(sympy.S(3)/4)/(15*b**2) - 2*x**3*(a - b*x**2)**(sympy.S(3)/4)/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_812():
    f = x**2/(a - b*x**2)**(sympy.S(1)/4)
    F = 4*a**(sympy.S(3)/2)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(5*b**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(1)/4)) - 2*x*(a - b*x**2)**(sympy.S(3)/4)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_813():
    f = (a - b*x**2)**(sympy.S(-1)/4)
    F = 2*sqrt(a)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(b)*(a - b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_814():
    f = 1/(x**2*(a - b*x**2)**(sympy.S(1)/4))
    F = -(a - b*x**2)**(sympy.S(3)/4)/(a*x) - sqrt(b)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(a)*(a - b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_815():
    f = 1/(x**4*(a - b*x**2)**(sympy.S(1)/4))
    F = -(a - b*x**2)**(sympy.S(3)/4)/(3*a*x**3) - b*(a - b*x**2)**(sympy.S(3)/4)/(2*a**2*x) - b**(sympy.S(3)/2)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(2*a**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_816():
    f = 1/(x**6*(a - b*x**2)**(sympy.S(1)/4))
    F = -(a - b*x**2)**(sympy.S(3)/4)/(5*a*x**5) - 7*b*(a - b*x**2)**(sympy.S(3)/4)/(30*a**2*x**3) - 7*b**2*(a - b*x**2)**(sympy.S(3)/4)/(20*a**3*x) - 7*b**(sympy.S(5)/2)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(20*a**(sympy.S(5)/2)*(a - b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_817():
    f = x**6/(a + b*x**2)**(sympy.S(3)/4)
    F = -80*a**(sympy.S(7)/2)*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(77*b**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(3)/4)) + 40*a**2*x*(a + b*x**2)**(sympy.S(1)/4)/(77*b**3) - 20*a*x**3*(a + b*x**2)**(sympy.S(1)/4)/(77*b**2) + 2*x**5*(a + b*x**2)**(sympy.S(1)/4)/(11*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_818():
    f = x**4/(a + b*x**2)**(sympy.S(3)/4)
    F = 8*a**(sympy.S(5)/2)*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(7*b**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/4)) - 4*a*x*(a + b*x**2)**(sympy.S(1)/4)/(7*b**2) + 2*x**3*(a + b*x**2)**(sympy.S(1)/4)/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_819():
    f = x**2/(a + b*x**2)**(sympy.S(3)/4)
    F = -4*a**(sympy.S(3)/2)*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(3*b**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/4)) + 2*x*(a + b*x**2)**(sympy.S(1)/4)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_820():
    f = (a + b*x**2)**(sympy.S(-3)/4)
    F = 2*sqrt(a)*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(b)*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_821():
    f = 1/(x**2*(a + b*x**2)**(sympy.S(3)/4))
    F = -(a + b*x**2)**(sympy.S(1)/4)/(a*x) - sqrt(b)*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(a)*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_822():
    f = 1/(x**4*(a + b*x**2)**(sympy.S(3)/4))
    F = -(a + b*x**2)**(sympy.S(1)/4)/(3*a*x**3) + 5*b*(a + b*x**2)**(sympy.S(1)/4)/(6*a**2*x) + 5*b**(sympy.S(3)/2)*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(6*a**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_823():
    f = 1/(x**6*(a + b*x**2)**(sympy.S(3)/4))
    F = -(a + b*x**2)**(sympy.S(1)/4)/(5*a*x**5) + 3*b*(a + b*x**2)**(sympy.S(1)/4)/(10*a**2*x**3) - 3*b**2*(a + b*x**2)**(sympy.S(1)/4)/(4*a**3*x) - 3*b**(sympy.S(5)/2)*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(4*a**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_824():
    f = x**6/(a - b*x**2)**(sympy.S(3)/4)
    F = 80*a**(sympy.S(7)/2)*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(77*b**(sympy.S(7)/2)*(a - b*x**2)**(sympy.S(3)/4)) - 40*a**2*x*(a - b*x**2)**(sympy.S(1)/4)/(77*b**3) - 20*a*x**3*(a - b*x**2)**(sympy.S(1)/4)/(77*b**2) - 2*x**5*(a - b*x**2)**(sympy.S(1)/4)/(11*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_825():
    f = x**4/(a - b*x**2)**(sympy.S(3)/4)
    F = 8*a**(sympy.S(5)/2)*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(7*b**(sympy.S(5)/2)*(a - b*x**2)**(sympy.S(3)/4)) - 4*a*x*(a - b*x**2)**(sympy.S(1)/4)/(7*b**2) - 2*x**3*(a - b*x**2)**(sympy.S(1)/4)/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_826():
    f = x**2/(a - b*x**2)**(sympy.S(3)/4)
    F = 4*a**(sympy.S(3)/2)*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(3*b**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(3)/4)) - 2*x*(a - b*x**2)**(sympy.S(1)/4)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_827():
    f = (a - b*x**2)**(sympy.S(-3)/4)
    F = 2*sqrt(a)*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(b)*(a - b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_828():
    f = 1/(x**2*(a - b*x**2)**(sympy.S(3)/4))
    F = -(a - b*x**2)**(sympy.S(1)/4)/(a*x) + sqrt(b)*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(a)*(a - b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_829():
    f = 1/(x**4*(a - b*x**2)**(sympy.S(3)/4))
    F = -(a - b*x**2)**(sympy.S(1)/4)/(3*a*x**3) - 5*b*(a - b*x**2)**(sympy.S(1)/4)/(6*a**2*x) + 5*b**(sympy.S(3)/2)*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(6*a**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_830():
    f = 1/(x**6*(a - b*x**2)**(sympy.S(3)/4))
    F = -(a - b*x**2)**(sympy.S(1)/4)/(5*a*x**5) - 3*b*(a - b*x**2)**(sympy.S(1)/4)/(10*a**2*x**3) - 3*b**2*(a - b*x**2)**(sympy.S(1)/4)/(4*a**3*x) + 3*b**(sympy.S(5)/2)*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(4*a**(sympy.S(5)/2)*(a - b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_831():
    f = x**6/(a + b*x**2)**(sympy.S(5)/4)
    F = -16*a**(sympy.S(5)/2)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(3*b**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(1)/4)) + 8*a**2*x/(3*b**3*(a + b*x**2)**(sympy.S(1)/4)) - 4*a*x**3/(9*b**2*(a + b*x**2)**(sympy.S(1)/4)) + 2*x**5/(9*b*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_832():
    f = x**4/(a + b*x**2)**(sympy.S(5)/4)
    F = 24*a**(sympy.S(3)/2)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(5*b**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)) - 12*a*x/(5*b**2*(a + b*x**2)**(sympy.S(1)/4)) + 2*x**3/(5*b*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_833():
    f = x**2/(a + b*x**2)**(sympy.S(5)/4)
    F = -4*sqrt(a)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(b**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)) + 2*x/(b*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_834():
    f = (a + b*x**2)**(sympy.S(-5)/4)
    F = 2*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(a)*sqrt(b)*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_835():
    f = 1/(x**2*(a + b*x**2)**(sympy.S(5)/4))
    F = -1/(a*x*(a + b*x**2)**(sympy.S(1)/4)) - 3*sqrt(b)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(a**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_836():
    f = 1/(x**4*(a + b*x**2)**(sympy.S(5)/4))
    F = -1/(3*a*x**3*(a + b*x**2)**(sympy.S(1)/4)) + 7*b/(6*a**2*x*(a + b*x**2)**(sympy.S(1)/4)) + 7*b**(sympy.S(3)/2)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(2*a**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_837():
    f = 1/(x**6*(a + b*x**2)**(sympy.S(5)/4))
    F = -1/(5*a*x**5*(a + b*x**2)**(sympy.S(1)/4)) + 11*b/(30*a**2*x**3*(a + b*x**2)**(sympy.S(1)/4)) - 77*b**2/(60*a**3*x*(a + b*x**2)**(sympy.S(1)/4)) - 77*b**(sympy.S(5)/2)*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(20*a**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_838():
    f = x**6/(a - b*x**2)**(sympy.S(5)/4)
    F = -16*a**(sympy.S(5)/2)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(3*b**(sympy.S(7)/2)*(a - b*x**2)**(sympy.S(1)/4)) + 8*a*x*(a - b*x**2)**(sympy.S(3)/4)/(3*b**3) + 2*x**5/(b*(a - b*x**2)**(sympy.S(1)/4)) + 20*x**3*(a - b*x**2)**(sympy.S(3)/4)/(9*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_839():
    f = x**4/(a - b*x**2)**(sympy.S(5)/4)
    F = -24*a**(sympy.S(3)/2)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(5*b**(sympy.S(5)/2)*(a - b*x**2)**(sympy.S(1)/4)) + 2*x**3/(b*(a - b*x**2)**(sympy.S(1)/4)) + 12*x*(a - b*x**2)**(sympy.S(3)/4)/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_840():
    f = x**2/(a - b*x**2)**(sympy.S(5)/4)
    F = -4*sqrt(a)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(b**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(1)/4)) + 2*x/(b*(a - b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_841():
    f = (a - b*x**2)**(sympy.S(-5)/4)
    F = 2*x/(a*(a - b*x**2)**(sympy.S(1)/4)) - 2*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(a)*sqrt(b)*(a - b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_842():
    f = 1/(x**2*(a - b*x**2)**(sympy.S(5)/4))
    F = 2/(a*x*(a - b*x**2)**(sympy.S(1)/4)) - 3*(a - b*x**2)**(sympy.S(3)/4)/(a**2*x) - 3*sqrt(b)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(a**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_843():
    f = 1/(x**4*(a - b*x**2)**(sympy.S(5)/4))
    F = 2/(a*x**3*(a - b*x**2)**(sympy.S(1)/4)) - 7*(a - b*x**2)**(sympy.S(3)/4)/(3*a**2*x**3) - 7*b*(a - b*x**2)**(sympy.S(3)/4)/(2*a**3*x) - 7*b**(sympy.S(3)/2)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(2*a**(sympy.S(5)/2)*(a - b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_844():
    f = 1/(x**6*(a - b*x**2)**(sympy.S(5)/4))
    F = 2/(a*x**5*(a - b*x**2)**(sympy.S(1)/4)) - 11*(a - b*x**2)**(sympy.S(3)/4)/(5*a**2*x**5) - 77*b*(a - b*x**2)**(sympy.S(3)/4)/(30*a**3*x**3) - 77*b**2*(a - b*x**2)**(sympy.S(3)/4)/(20*a**4*x) - 77*b**(sympy.S(5)/2)*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(20*a**(sympy.S(7)/2)*(a - b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_845():
    f = (a + b*x**2)**(sympy.S(-7)/4)
    F = 2*x/(3*a*(a + b*x**2)**(sympy.S(3)/4)) + 2*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(3*sqrt(a)*sqrt(b)*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_846():
    f = (a + b*x**2)**(sympy.S(-9)/4)
    F = 2*x/(5*a*(a + b*x**2)**(sympy.S(5)/4)) + 6*(1 + b*x**2/a)**(sympy.S(1)/4)*elliptic_e(atan(sqrt(b)*x/sqrt(a))/2, 2)/(5*a**(sympy.S(3)/2)*sqrt(b)*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_847():
    f = (a + b*x**2)**(sympy.S(-11)/4)
    F = 2*x/(7*a*(a + b*x**2)**(sympy.S(7)/4)) + 10*x/(21*a**2*(a + b*x**2)**(sympy.S(3)/4)) + 10*(1 + b*x**2/a)**(sympy.S(3)/4)*elliptic_f(atan(sqrt(b)*x/sqrt(a))/2, 2)/(21*a**(sympy.S(3)/2)*sqrt(b)*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_848():
    f = (a - b*x**2)**(sympy.S(-7)/4)
    F = 2*x/(3*a*(a - b*x**2)**(sympy.S(3)/4)) + 2*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(3*sqrt(a)*sqrt(b)*(a - b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_849():
    f = (a - b*x**2)**(sympy.S(-9)/4)
    F = 2*x/(5*a*(a - b*x**2)**(sympy.S(5)/4)) + 6*x/(5*a**2*(a - b*x**2)**(sympy.S(1)/4)) - 6*(1 - b*x**2/a)**(sympy.S(1)/4)*elliptic_e(asin(sqrt(b)*x/sqrt(a))/2, 2)/(5*a**(sympy.S(3)/2)*sqrt(b)*(a - b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_850():
    f = (a - b*x**2)**(sympy.S(-11)/4)
    F = 2*x/(7*a*(a - b*x**2)**(sympy.S(7)/4)) + 10*x/(21*a**2*(a - b*x**2)**(sympy.S(3)/4)) + 10*(1 - b*x**2/a)**(sympy.S(3)/4)*elliptic_f(asin(sqrt(b)*x/sqrt(a))/2, 2)/(21*a**(sympy.S(3)/2)*sqrt(b)*(a - b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_851():
    f = x**6/(3*x**2 + 2)**(sympy.S(1)/4)
    F = 2*x**5*(3*x**2 + 2)**(sympy.S(3)/4)/39 - 40*x**3*(3*x**2 + 2)**(sympy.S(3)/4)/1053 + 32*x*(3*x**2 + 2)**(sympy.S(3)/4)/1053 - 128*x/(1053*(3*x**2 + 2)**(sympy.S(1)/4)) + 128*2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(atan(sqrt(6)*x/2)/2, 2)/3159
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_852():
    f = x**4/(3*x**2 + 2)**(sympy.S(1)/4)
    F = 2*x**3*(3*x**2 + 2)**(sympy.S(3)/4)/27 - 8*x*(3*x**2 + 2)**(sympy.S(3)/4)/135 + 32*x/(135*(3*x**2 + 2)**(sympy.S(1)/4)) - 32*2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(atan(sqrt(6)*x/2)/2, 2)/405
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_853():
    f = x**2/(3*x**2 + 2)**(sympy.S(1)/4)
    F = 2*x*(3*x**2 + 2)**(sympy.S(3)/4)/15 - 8*x/(15*(3*x**2 + 2)**(sympy.S(1)/4)) + 8*2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(atan(sqrt(6)*x/2)/2, 2)/45
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_854():
    f = (3*x**2 + 2)**(sympy.S(-1)/4)
    F = 2*x/(3*x**2 + 2)**(sympy.S(1)/4) - 2*2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(atan(sqrt(6)*x/2)/2, 2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_855():
    f = 1/(x**2*(3*x**2 + 2)**(sympy.S(1)/4))
    F = 3*x/(2*(3*x**2 + 2)**(sympy.S(1)/4)) - 2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(atan(sqrt(6)*x/2)/2, 2)/2 - (3*x**2 + 2)**(sympy.S(3)/4)/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_856():
    f = 1/(x**4*(3*x**2 + 2)**(sympy.S(1)/4))
    F = -9*x/(8*(3*x**2 + 2)**(sympy.S(1)/4)) + 3*2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(atan(sqrt(6)*x/2)/2, 2)/8 + 3*(3*x**2 + 2)**(sympy.S(3)/4)/(8*x) - (3*x**2 + 2)**(sympy.S(3)/4)/(6*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_857():
    f = 1/(x**6*(3*x**2 + 2)**(sympy.S(1)/4))
    F = 189*x/(160*(3*x**2 + 2)**(sympy.S(1)/4)) - 63*2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(atan(sqrt(6)*x/2)/2, 2)/160 - 63*(3*x**2 + 2)**(sympy.S(3)/4)/(160*x) + 7*(3*x**2 + 2)**(sympy.S(3)/4)/(40*x**3) - (3*x**2 + 2)**(sympy.S(3)/4)/(10*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_858():
    f = x**6/(2 - 3*x**2)**(sympy.S(1)/4)
    F = -2*x**5*(2 - 3*x**2)**(sympy.S(3)/4)/39 - 40*x**3*(2 - 3*x**2)**(sympy.S(3)/4)/1053 - 32*x*(2 - 3*x**2)**(sympy.S(3)/4)/1053 + 128*2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(asin(sqrt(6)*x/2)/2, 2)/3159
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_859():
    f = x**4/(2 - 3*x**2)**(sympy.S(1)/4)
    F = -2*x**3*(2 - 3*x**2)**(sympy.S(3)/4)/27 - 8*x*(2 - 3*x**2)**(sympy.S(3)/4)/135 + 32*2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(asin(sqrt(6)*x/2)/2, 2)/405
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_860():
    f = x**2/(2 - 3*x**2)**(sympy.S(1)/4)
    F = -2*x*(2 - 3*x**2)**(sympy.S(3)/4)/15 + 8*2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(asin(sqrt(6)*x/2)/2, 2)/45
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_861():
    f = (2 - 3*x**2)**(sympy.S(-1)/4)
    F = 2*2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(asin(sqrt(6)*x/2)/2, 2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_862():
    f = 1/(x**2*(2 - 3*x**2)**(sympy.S(1)/4))
    F = -2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(asin(sqrt(6)*x/2)/2, 2)/2 - (2 - 3*x**2)**(sympy.S(3)/4)/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_863():
    f = 1/(x**4*(2 - 3*x**2)**(sympy.S(1)/4))
    F = -3*2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(asin(sqrt(6)*x/2)/2, 2)/8 - 3*(2 - 3*x**2)**(sympy.S(3)/4)/(8*x) - (2 - 3*x**2)**(sympy.S(3)/4)/(6*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_864():
    f = 1/(x**6*(2 - 3*x**2)**(sympy.S(1)/4))
    F = -63*2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(asin(sqrt(6)*x/2)/2, 2)/160 - 63*(2 - 3*x**2)**(sympy.S(3)/4)/(160*x) - 7*(2 - 3*x**2)**(sympy.S(3)/4)/(40*x**3) - (2 - 3*x**2)**(sympy.S(3)/4)/(10*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_865():
    f = x**6/(3*x**2 + 2)**(sympy.S(3)/4)
    F = 2*x**5*(3*x**2 + 2)**(sympy.S(1)/4)/33 - 40*x**3*(3*x**2 + 2)**(sympy.S(1)/4)/693 + 160*x*(3*x**2 + 2)**(sympy.S(1)/4)/2079 - 320*2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(atan(sqrt(6)*x/2)/2, 2)/6237
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_866():
    f = x**4/(3*x**2 + 2)**(sympy.S(3)/4)
    F = 2*x**3*(3*x**2 + 2)**(sympy.S(1)/4)/21 - 8*x*(3*x**2 + 2)**(sympy.S(1)/4)/63 + 16*2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(atan(sqrt(6)*x/2)/2, 2)/189
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_867():
    f = x**2/(3*x**2 + 2)**(sympy.S(3)/4)
    F = 2*x*(3*x**2 + 2)**(sympy.S(1)/4)/9 - 4*2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(atan(sqrt(6)*x/2)/2, 2)/27
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_868():
    f = (3*x**2 + 2)**(sympy.S(-3)/4)
    F = 2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(atan(sqrt(6)*x/2)/2, 2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_869():
    f = 1/(x**2*(3*x**2 + 2)**(sympy.S(3)/4))
    F = -2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(atan(sqrt(6)*x/2)/2, 2)/4 - (3*x**2 + 2)**(sympy.S(1)/4)/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_870():
    f = 1/(x**4*(3*x**2 + 2)**(sympy.S(3)/4))
    F = 5*2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(atan(sqrt(6)*x/2)/2, 2)/16 + 5*(3*x**2 + 2)**(sympy.S(1)/4)/(8*x) - (3*x**2 + 2)**(sympy.S(1)/4)/(6*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_871():
    f = 1/(x**6*(3*x**2 + 2)**(sympy.S(3)/4))
    F = -27*2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(atan(sqrt(6)*x/2)/2, 2)/64 - 27*(3*x**2 + 2)**(sympy.S(1)/4)/(32*x) + 9*(3*x**2 + 2)**(sympy.S(1)/4)/(40*x**3) - (3*x**2 + 2)**(sympy.S(1)/4)/(10*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_872():
    f = x**6/(2 - 3*x**2)**(sympy.S(3)/4)
    F = -2*x**5*(2 - 3*x**2)**(sympy.S(1)/4)/33 - 40*x**3*(2 - 3*x**2)**(sympy.S(1)/4)/693 - 160*x*(2 - 3*x**2)**(sympy.S(1)/4)/2079 + 320*2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(asin(sqrt(6)*x/2)/2, 2)/6237
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_873():
    f = x**4/(2 - 3*x**2)**(sympy.S(3)/4)
    F = -2*x**3*(2 - 3*x**2)**(sympy.S(1)/4)/21 - 8*x*(2 - 3*x**2)**(sympy.S(1)/4)/63 + 16*2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(asin(sqrt(6)*x/2)/2, 2)/189
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_874():
    f = x**2/(2 - 3*x**2)**(sympy.S(3)/4)
    F = -2*x*(2 - 3*x**2)**(sympy.S(1)/4)/9 + 4*2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(asin(sqrt(6)*x/2)/2, 2)/27
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_875():
    f = (2 - 3*x**2)**(sympy.S(-3)/4)
    F = 2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(asin(sqrt(6)*x/2)/2, 2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_876():
    f = 1/(x**2*(2 - 3*x**2)**(sympy.S(3)/4))
    F = 2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(asin(sqrt(6)*x/2)/2, 2)/4 - (2 - 3*x**2)**(sympy.S(1)/4)/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_877():
    f = 1/(x**4*(2 - 3*x**2)**(sympy.S(3)/4))
    F = 5*2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(asin(sqrt(6)*x/2)/2, 2)/16 - 5*(2 - 3*x**2)**(sympy.S(1)/4)/(8*x) - (2 - 3*x**2)**(sympy.S(1)/4)/(6*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_878():
    f = 1/(x**6*(2 - 3*x**2)**(sympy.S(3)/4))
    F = 27*2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(asin(sqrt(6)*x/2)/2, 2)/64 - 27*(2 - 3*x**2)**(sympy.S(1)/4)/(32*x) - 9*(2 - 3*x**2)**(sympy.S(1)/4)/(40*x**3) - (2 - 3*x**2)**(sympy.S(1)/4)/(10*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_879():
    f = x**6/(3*x**2 - 2)**(sympy.S(1)/4)
    F = 2*x**5*(3*x**2 - 2)**(sympy.S(3)/4)/39 + 40*x**3*(3*x**2 - 2)**(sympy.S(3)/4)/1053 + 32*x*(3*x**2 - 2)**(sympy.S(3)/4)/1053 + 128*x*(3*x**2 - 2)**(sympy.S(1)/4)/(1053*sqrt(3*x**2 - 2) + 1053*sqrt(2)) - 128*2**(sympy.S(1)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_e(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(3159*x) + 64*2**(sympy.S(1)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(3159*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_880():
    f = x**4/(3*x**2 - 2)**(sympy.S(1)/4)
    F = 2*x**3*(3*x**2 - 2)**(sympy.S(3)/4)/27 + 8*x*(3*x**2 - 2)**(sympy.S(3)/4)/135 + 32*x*(3*x**2 - 2)**(sympy.S(1)/4)/(135*sqrt(3*x**2 - 2) + 135*sqrt(2)) - 32*2**(sympy.S(1)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_e(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(405*x) + 16*2**(sympy.S(1)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(405*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_881():
    f = x**2/(3*x**2 - 2)**(sympy.S(1)/4)
    F = 2*x*(3*x**2 - 2)**(sympy.S(3)/4)/15 + 8*x*(3*x**2 - 2)**(sympy.S(1)/4)/(15*sqrt(3*x**2 - 2) + 15*sqrt(2)) - 8*2**(sympy.S(1)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_e(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(45*x) + 4*2**(sympy.S(1)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(45*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_882():
    f = (3*x**2 - 2)**(sympy.S(-1)/4)
    F = 2*x*(3*x**2 - 2)**(sympy.S(1)/4)/(sqrt(3*x**2 - 2) + sqrt(2)) - 2*2**(sympy.S(1)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_e(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(3*x) + 2**(sympy.S(1)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_883():
    f = 1/(x**2*(3*x**2 - 2)**(sympy.S(1)/4))
    F = -3*x*(3*x**2 - 2)**(sympy.S(1)/4)/(2*sqrt(3*x**2 - 2) + 2*sqrt(2)) + 2**(sympy.S(1)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_e(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(2*x) - 2**(sympy.S(1)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(4*x) + (3*x**2 - 2)**(sympy.S(3)/4)/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_884():
    f = 1/(x**4*(3*x**2 - 2)**(sympy.S(1)/4))
    F = -9*x*(3*x**2 - 2)**(sympy.S(1)/4)/(8*sqrt(3*x**2 - 2) + 8*sqrt(2)) + 3*2**(sympy.S(1)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_e(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(8*x) - 3*2**(sympy.S(1)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(16*x) + 3*(3*x**2 - 2)**(sympy.S(3)/4)/(8*x) + (3*x**2 - 2)**(sympy.S(3)/4)/(6*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_885():
    f = 1/(x**6*(3*x**2 - 2)**(sympy.S(1)/4))
    F = -189*x*(3*x**2 - 2)**(sympy.S(1)/4)/(160*sqrt(3*x**2 - 2) + 160*sqrt(2)) + 63*2**(sympy.S(1)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_e(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(160*x) - 63*2**(sympy.S(1)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(320*x) + 63*(3*x**2 - 2)**(sympy.S(3)/4)/(160*x) + 7*(3*x**2 - 2)**(sympy.S(3)/4)/(40*x**3) + (3*x**2 - 2)**(sympy.S(3)/4)/(10*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_886():
    f = x**6/(-3*x**2 - 2)**(sympy.S(1)/4)
    F = -2*x**5*(-3*x**2 - 2)**(sympy.S(3)/4)/39 + 40*x**3*(-3*x**2 - 2)**(sympy.S(3)/4)/1053 - 32*x*(-3*x**2 - 2)**(sympy.S(3)/4)/1053 - 128*x*(-3*x**2 - 2)**(sympy.S(1)/4)/(1053*sqrt(-3*x**2 - 2) + 1053*sqrt(2)) - 128*2**(sympy.S(1)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_e(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(3159*x) + 64*2**(sympy.S(1)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(3159*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_887():
    f = x**4/(-3*x**2 - 2)**(sympy.S(1)/4)
    F = -2*x**3*(-3*x**2 - 2)**(sympy.S(3)/4)/27 + 8*x*(-3*x**2 - 2)**(sympy.S(3)/4)/135 + 32*x*(-3*x**2 - 2)**(sympy.S(1)/4)/(135*sqrt(-3*x**2 - 2) + 135*sqrt(2)) + 32*2**(sympy.S(1)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_e(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(405*x) - 16*2**(sympy.S(1)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(405*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_888():
    f = x**2/(-3*x**2 - 2)**(sympy.S(1)/4)
    F = -2*x*(-3*x**2 - 2)**(sympy.S(3)/4)/15 - 8*x*(-3*x**2 - 2)**(sympy.S(1)/4)/(15*sqrt(-3*x**2 - 2) + 15*sqrt(2)) - 8*2**(sympy.S(1)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_e(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(45*x) + 4*2**(sympy.S(1)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(45*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_889():
    f = (-3*x**2 - 2)**(sympy.S(-1)/4)
    F = 2*x*(-3*x**2 - 2)**(sympy.S(1)/4)/(sqrt(-3*x**2 - 2) + sqrt(2)) + 2*2**(sympy.S(1)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_e(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(3*x) - 2**(sympy.S(1)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_890():
    f = 1/(x**2*(-3*x**2 - 2)**(sympy.S(1)/4))
    F = 3*x*(-3*x**2 - 2)**(sympy.S(1)/4)/(2*sqrt(-3*x**2 - 2) + 2*sqrt(2)) + 2**(sympy.S(1)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_e(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(2*x) - 2**(sympy.S(1)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(4*x) + (-3*x**2 - 2)**(sympy.S(3)/4)/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_891():
    f = 1/(x**4*(-3*x**2 - 2)**(sympy.S(1)/4))
    F = -9*x*(-3*x**2 - 2)**(sympy.S(1)/4)/(8*sqrt(-3*x**2 - 2) + 8*sqrt(2)) - 3*2**(sympy.S(1)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_e(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(8*x) + 3*2**(sympy.S(1)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(16*x) - 3*(-3*x**2 - 2)**(sympy.S(3)/4)/(8*x) + (-3*x**2 - 2)**(sympy.S(3)/4)/(6*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_892():
    f = 1/(x**6*(-3*x**2 - 2)**(sympy.S(1)/4))
    F = 189*x*(-3*x**2 - 2)**(sympy.S(1)/4)/(160*sqrt(-3*x**2 - 2) + 160*sqrt(2)) + 63*2**(sympy.S(1)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_e(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(160*x) - 63*2**(sympy.S(1)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(320*x) + 63*(-3*x**2 - 2)**(sympy.S(3)/4)/(160*x) - 7*(-3*x**2 - 2)**(sympy.S(3)/4)/(40*x**3) + (-3*x**2 - 2)**(sympy.S(3)/4)/(10*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_893():
    f = x**6/(3*x**2 - 2)**(sympy.S(3)/4)
    F = 2*x**5*(3*x**2 - 2)**(sympy.S(1)/4)/33 + 40*x**3*(3*x**2 - 2)**(sympy.S(1)/4)/693 + 160*x*(3*x**2 - 2)**(sympy.S(1)/4)/2079 + 160*2**(sympy.S(3)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(6237*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_894():
    f = x**4/(3*x**2 - 2)**(sympy.S(3)/4)
    F = 2*x**3*(3*x**2 - 2)**(sympy.S(1)/4)/21 + 8*x*(3*x**2 - 2)**(sympy.S(1)/4)/63 + 8*2**(sympy.S(3)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(189*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_895():
    f = x**2/(3*x**2 - 2)**(sympy.S(3)/4)
    F = 2*x*(3*x**2 - 2)**(sympy.S(1)/4)/9 + 2*2**(sympy.S(3)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(27*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_896():
    f = (3*x**2 - 2)**(sympy.S(-3)/4)
    F = 2**(sympy.S(3)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(6*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_897():
    f = 1/(x**2*(3*x**2 - 2)**(sympy.S(3)/4))
    F = 2**(sympy.S(3)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(8*x) + (3*x**2 - 2)**(sympy.S(1)/4)/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_898():
    f = 1/(x**4*(3*x**2 - 2)**(sympy.S(3)/4))
    F = 5*2**(sympy.S(3)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(32*x) + 5*(3*x**2 - 2)**(sympy.S(1)/4)/(8*x) + (3*x**2 - 2)**(sympy.S(1)/4)/(6*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_899():
    f = 1/(x**6*(3*x**2 - 2)**(sympy.S(3)/4))
    F = 27*2**(sympy.S(3)/4)*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 2) + sqrt(2))**2)*(sqrt(3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(128*x) + 27*(3*x**2 - 2)**(sympy.S(1)/4)/(32*x) + 9*(3*x**2 - 2)**(sympy.S(1)/4)/(40*x**3) + (3*x**2 - 2)**(sympy.S(1)/4)/(10*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_900():
    f = x**6/(-3*x**2 - 2)**(sympy.S(3)/4)
    F = -2*x**5*(-3*x**2 - 2)**(sympy.S(1)/4)/33 + 40*x**3*(-3*x**2 - 2)**(sympy.S(1)/4)/693 - 160*x*(-3*x**2 - 2)**(sympy.S(1)/4)/2079 + 160*2**(sympy.S(3)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(6237*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_901():
    f = x**4/(-3*x**2 - 2)**(sympy.S(3)/4)
    F = -2*x**3*(-3*x**2 - 2)**(sympy.S(1)/4)/21 + 8*x*(-3*x**2 - 2)**(sympy.S(1)/4)/63 - 8*2**(sympy.S(3)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(189*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_902():
    f = x**2/(-3*x**2 - 2)**(sympy.S(3)/4)
    F = -2*x*(-3*x**2 - 2)**(sympy.S(1)/4)/9 + 2*2**(sympy.S(3)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(27*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_903():
    f = (-3*x**2 - 2)**(sympy.S(-3)/4)
    F = -2**(sympy.S(3)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(6*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_904():
    f = 1/(x**2*(-3*x**2 - 2)**(sympy.S(3)/4))
    F = 2**(sympy.S(3)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(8*x) + (-3*x**2 - 2)**(sympy.S(1)/4)/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_905():
    f = 1/(x**4*(-3*x**2 - 2)**(sympy.S(3)/4))
    F = -5*2**(sympy.S(3)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(32*x) - 5*(-3*x**2 - 2)**(sympy.S(1)/4)/(8*x) + (-3*x**2 - 2)**(sympy.S(1)/4)/(6*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_906():
    f = 1/(x**6*(-3*x**2 - 2)**(sympy.S(3)/4))
    F = 27*2**(sympy.S(3)/4)*sqrt(3)*sqrt(-x**2/(sqrt(-3*x**2 - 2) + sqrt(2))**2)*(sqrt(-3*x**2 - 2) + sqrt(2))*elliptic_f(2*atan(2**(sympy.S(3)/4)*(-3*x**2 - 2)**(sympy.S(1)/4)/2), sympy.S.Half)/(128*x) + 27*(-3*x**2 - 2)**(sympy.S(1)/4)/(32*x) - 9*(-3*x**2 - 2)**(sympy.S(1)/4)/(40*x**3) + (-3*x**2 - 2)**(sympy.S(1)/4)/(10*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_907():
    f = (c*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(1)/4)
    F = -a**(sympy.S(5)/2)*c**2*(c*x)**(sympy.S(3)/2)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(12*b**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/4)) - a**2*c**3*sqrt(c*x)*(a + b*x**2)**(sympy.S(1)/4)/(12*b**2) + a*c*(c*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)/(30*b) + (c*x)**(sympy.S(9)/2)*(a + b*x**2)**(sympy.S(1)/4)/(5*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_908():
    f = (c*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)
    F = a**(sympy.S(3)/2)*(c*x)**(sympy.S(3)/2)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(6*sqrt(b)*(a + b*x**2)**(sympy.S(3)/4)) + a*c*sqrt(c*x)*(a + b*x**2)**(sympy.S(1)/4)/(6*b) + (c*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_909():
    f = (a + b*x**2)**(sympy.S(1)/4)/sqrt(c*x)
    F = -sqrt(a)*sqrt(b)*(c*x)**(sympy.S(3)/2)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(c**2*(a + b*x**2)**(sympy.S(3)/4)) + sqrt(c*x)*(a + b*x**2)**(sympy.S(1)/4)/c
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_910():
    f = (a + b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(5)/2)
    F = -2*(a + b*x**2)**(sympy.S(1)/4)/(3*c*(c*x)**(sympy.S(3)/2)) - 2*b**(sympy.S(3)/2)*(c*x)**(sympy.S(3)/2)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(3*sqrt(a)*c**4*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_911():
    f = (a + b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(9)/2)
    F = -2*(a + b*x**2)**(sympy.S(1)/4)/(7*c*(c*x)**(sympy.S(7)/2)) - 2*b*(a + b*x**2)**(sympy.S(1)/4)/(21*a*c**3*(c*x)**(sympy.S(3)/2)) + 4*b**(sympy.S(5)/2)*(c*x)**(sympy.S(3)/2)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(21*a**(sympy.S(3)/2)*c**6*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_912():
    f = (a + b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(13)/2)
    F = -2*(a + b*x**2)**(sympy.S(1)/4)/(11*c*(c*x)**(sympy.S(11)/2)) - 2*b*(a + b*x**2)**(sympy.S(1)/4)/(77*a*c**3*(c*x)**(sympy.S(7)/2)) + 4*b**2*(a + b*x**2)**(sympy.S(1)/4)/(77*a**2*c**5*(c*x)**(sympy.S(3)/2)) - 8*b**(sympy.S(7)/2)*(c*x)**(sympy.S(3)/2)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(77*a**(sympy.S(5)/2)*c**8*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_913():
    f = (c*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)
    F = 3*a**2*c**(sympy.S(5)/2)*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/(32*b**(sympy.S(7)/4)) - 3*a**2*c**(sympy.S(5)/2)*atanh(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/(32*b**(sympy.S(7)/4)) + a*c*(c*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)/(16*b) + (c*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(1)/4)/(4*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_914():
    f = sqrt(c*x)*(a + b*x**2)**(sympy.S(1)/4)
    F = -a*sqrt(c)*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(3)/4)) + a*sqrt(c)*atanh(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(3)/4)) + (c*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)/(2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_915():
    f = (a + b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(3)/2)
    F = -b**(sympy.S(1)/4)*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/c**(sympy.S(3)/2) + b**(sympy.S(1)/4)*atanh(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/c**(sympy.S(3)/2) - 2*(a + b*x**2)**(sympy.S(1)/4)/(c*sqrt(c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_916():
    f = (a + b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(7)/2)
    F = -2*(a + b*x**2)**(sympy.S(5)/4)/(5*a*c*(c*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_917():
    f = (a + b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(11)/2)
    F = -2*(a + b*x**2)**(sympy.S(5)/4)/(5*a*c*(c*x)**(sympy.S(9)/2)) + 8*(a + b*x**2)**(sympy.S(9)/4)/(45*a**2*c*(c*x)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_918():
    f = (a + b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(15)/2)
    F = -2*(a + b*x**2)**(sympy.S(5)/4)/(5*a*c*(c*x)**(sympy.S(13)/2)) + 16*(a + b*x**2)**(sympy.S(9)/4)/(45*a**2*c*(c*x)**(sympy.S(13)/2)) - 64*(a + b*x**2)**(sympy.S(13)/4)/(585*a**3*c*(c*x)**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_919():
    f = (a + b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(19)/2)
    F = -2*(a + b*x**2)**(sympy.S(5)/4)/(5*a*c*(c*x)**(sympy.S(17)/2)) + 8*(a + b*x**2)**(sympy.S(9)/4)/(15*a**2*c*(c*x)**(sympy.S(17)/2)) - 64*(a + b*x**2)**(sympy.S(13)/4)/(195*a**3*c*(c*x)**(sympy.S(17)/2)) + 256*(a + b*x**2)**(sympy.S(17)/4)/(3315*a**4*c*(c*x)**(sympy.S(17)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_920():
    f = (c*x)**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(1)/4)
    F = -a**(sympy.S(3)/2)*(c*x)**(sympy.S(3)/2)*(-a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acsc(sqrt(b)*x/sqrt(a))/2, 2)/(6*sqrt(b)*(a - b*x**2)**(sympy.S(3)/4)) - a*c*sqrt(c*x)*(a - b*x**2)**(sympy.S(1)/4)/(6*b) + (c*x)**(sympy.S(5)/2)*(a - b*x**2)**(sympy.S(1)/4)/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_921():
    f = (a - b*x**2)**(sympy.S(1)/4)/sqrt(c*x)
    F = -sqrt(a)*sqrt(b)*(c*x)**(sympy.S(3)/2)*(-a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acsc(sqrt(b)*x/sqrt(a))/2, 2)/(c**2*(a - b*x**2)**(sympy.S(3)/4)) + sqrt(c*x)*(a - b*x**2)**(sympy.S(1)/4)/c
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_922():
    f = (a - b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(5)/2)
    F = -2*(a - b*x**2)**(sympy.S(1)/4)/(3*c*(c*x)**(sympy.S(3)/2)) + 2*b**(sympy.S(3)/2)*(c*x)**(sympy.S(3)/2)*(-a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acsc(sqrt(b)*x/sqrt(a))/2, 2)/(3*sqrt(a)*c**4*(a - b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_923():
    f = (a - b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(9)/2)
    F = -2*(a - b*x**2)**(sympy.S(1)/4)/(7*c*(c*x)**(sympy.S(7)/2)) + 2*b*(a - b*x**2)**(sympy.S(1)/4)/(21*a*c**3*(c*x)**(sympy.S(3)/2)) + 4*b**(sympy.S(5)/2)*(c*x)**(sympy.S(3)/2)*(-a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acsc(sqrt(b)*x/sqrt(a))/2, 2)/(21*a**(sympy.S(3)/2)*c**6*(a - b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_924():
    f = (a - b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(13)/2)
    F = -2*(a - b*x**2)**(sympy.S(1)/4)/(11*c*(c*x)**(sympy.S(11)/2)) + 2*b*(a - b*x**2)**(sympy.S(1)/4)/(77*a*c**3*(c*x)**(sympy.S(7)/2)) + 4*b**2*(a - b*x**2)**(sympy.S(1)/4)/(77*a**2*c**5*(c*x)**(sympy.S(3)/2)) + 8*b**(sympy.S(7)/2)*(c*x)**(sympy.S(3)/2)*(-a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acsc(sqrt(b)*x/sqrt(a))/2, 2)/(77*a**(sympy.S(5)/2)*c**8*(a - b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_925():
    f = (c*x)**(sympy.S(5)/2)*(a - b*x**2)**(sympy.S(1)/4)
    F = 3*sqrt(2)*a**2*c**(sympy.S(5)/2)*log(-sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(a - b*x**2)**(sympy.S(1)/4) + sqrt(b)*sqrt(c)*x/sqrt(a - b*x**2) + sqrt(c))/(128*b**(sympy.S(7)/4)) - 3*sqrt(2)*a**2*c**(sympy.S(5)/2)*log(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(a - b*x**2)**(sympy.S(1)/4) + sqrt(b)*sqrt(c)*x/sqrt(a - b*x**2) + sqrt(c))/(128*b**(sympy.S(7)/4)) + 3*sqrt(2)*a**2*c**(sympy.S(5)/2)*atan(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a - b*x**2)**(sympy.S(1)/4)) - 1)/(64*b**(sympy.S(7)/4)) + 3*sqrt(2)*a**2*c**(sympy.S(5)/2)*atan(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a - b*x**2)**(sympy.S(1)/4)) + 1)/(64*b**(sympy.S(7)/4)) - a*c*(c*x)**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(1)/4)/(16*b) + (c*x)**(sympy.S(7)/2)*(a - b*x**2)**(sympy.S(1)/4)/(4*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_926():
    f = sqrt(c*x)*(a - b*x**2)**(sympy.S(1)/4)
    F = sqrt(2)*a*sqrt(c)*log(-sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(a - b*x**2)**(sympy.S(1)/4) + sqrt(b)*sqrt(c)*x/sqrt(a - b*x**2) + sqrt(c))/(16*b**(sympy.S(3)/4)) - sqrt(2)*a*sqrt(c)*log(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(a - b*x**2)**(sympy.S(1)/4) + sqrt(b)*sqrt(c)*x/sqrt(a - b*x**2) + sqrt(c))/(16*b**(sympy.S(3)/4)) + sqrt(2)*a*sqrt(c)*atan(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a - b*x**2)**(sympy.S(1)/4)) - 1)/(8*b**(sympy.S(3)/4)) + sqrt(2)*a*sqrt(c)*atan(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a - b*x**2)**(sympy.S(1)/4)) + 1)/(8*b**(sympy.S(3)/4)) + (c*x)**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(1)/4)/(2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_927():
    f = (a - b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(3)/2)
    F = -sqrt(2)*b**(sympy.S(1)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(a - b*x**2)**(sympy.S(1)/4) + sqrt(b)*sqrt(c)*x/sqrt(a - b*x**2) + sqrt(c))/(4*c**(sympy.S(3)/2)) + sqrt(2)*b**(sympy.S(1)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(a - b*x**2)**(sympy.S(1)/4) + sqrt(b)*sqrt(c)*x/sqrt(a - b*x**2) + sqrt(c))/(4*c**(sympy.S(3)/2)) - sqrt(2)*b**(sympy.S(1)/4)*atan(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a - b*x**2)**(sympy.S(1)/4)) - 1)/(2*c**(sympy.S(3)/2)) - sqrt(2)*b**(sympy.S(1)/4)*atan(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a - b*x**2)**(sympy.S(1)/4)) + 1)/(2*c**(sympy.S(3)/2)) - 2*(a - b*x**2)**(sympy.S(1)/4)/(c*sqrt(c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_928():
    f = (a - b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(7)/2)
    F = -2*(a - b*x**2)**(sympy.S(5)/4)/(5*a*c*(c*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_929():
    f = (a - b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(11)/2)
    F = -2*(a - b*x**2)**(sympy.S(5)/4)/(5*a*c*(c*x)**(sympy.S(9)/2)) + 8*(a - b*x**2)**(sympy.S(9)/4)/(45*a**2*c*(c*x)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_930():
    f = (a - b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(15)/2)
    F = -2*(a - b*x**2)**(sympy.S(5)/4)/(5*a*c*(c*x)**(sympy.S(13)/2)) + 16*(a - b*x**2)**(sympy.S(9)/4)/(45*a**2*c*(c*x)**(sympy.S(13)/2)) - 64*(a - b*x**2)**(sympy.S(13)/4)/(585*a**3*c*(c*x)**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_931():
    f = (a - b*x**2)**(sympy.S(1)/4)/(c*x)**(sympy.S(19)/2)
    F = -2*(a - b*x**2)**(sympy.S(5)/4)/(5*a*c*(c*x)**(sympy.S(17)/2)) + 8*(a - b*x**2)**(sympy.S(9)/4)/(15*a**2*c*(c*x)**(sympy.S(17)/2)) - 64*(a - b*x**2)**(sympy.S(13)/4)/(195*a**3*c*(c*x)**(sympy.S(17)/2)) + 256*(a - b*x**2)**(sympy.S(17)/4)/(3315*a**4*c*(c*x)**(sympy.S(17)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_932():
    f = (c*x)**(sympy.S(3)/2)/(a + b*x**2)**(sympy.S(1)/4)
    F = -a*c**(sympy.S(3)/2)*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(5)/4)) - a*c**(sympy.S(3)/2)*atanh(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(5)/4)) + c*sqrt(c*x)*(a + b*x**2)**(sympy.S(3)/4)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_933():
    f = 1/(sqrt(c*x)*(a + b*x**2)**(sympy.S(1)/4))
    F = atan(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/(b**(sympy.S(1)/4)*sqrt(c)) + atanh(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/(b**(sympy.S(1)/4)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_934():
    f = 1/((c*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4))
    F = -2*(a + b*x**2)**(sympy.S(3)/4)/(3*a*c*(c*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_935():
    f = 1/((c*x)**(sympy.S(9)/2)*(a + b*x**2)**(sympy.S(1)/4))
    F = -2*(a + b*x**2)**(sympy.S(3)/4)/(3*a*c*(c*x)**(sympy.S(7)/2)) + 8*(a + b*x**2)**(sympy.S(7)/4)/(21*a**2*c*(c*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_936():
    f = 1/((c*x)**(sympy.S(13)/2)*(a + b*x**2)**(sympy.S(1)/4))
    F = -2*(a + b*x**2)**(sympy.S(3)/4)/(3*a*c*(c*x)**(sympy.S(11)/2)) + 16*(a + b*x**2)**(sympy.S(7)/4)/(21*a**2*c*(c*x)**(sympy.S(11)/2)) - 64*(a + b*x**2)**(sympy.S(11)/4)/(231*a**3*c*(c*x)**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_937():
    f = (c*x)**(sympy.S(9)/2)/(a + b*x**2)**(sympy.S(1)/4)
    F = 7*a**(sympy.S(5)/2)*c**4*sqrt(c*x)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(20*b**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)) + 7*a**2*c**4*x*sqrt(c*x)/(20*b**2*(a + b*x**2)**(sympy.S(1)/4)) - 7*a*c**3*(c*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/4)/(30*b**2) + c*(c*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(3)/4)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_938():
    f = (c*x)**(sympy.S(5)/2)/(a + b*x**2)**(sympy.S(1)/4)
    F = -a**(sympy.S(3)/2)*c**2*sqrt(c*x)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(2*b**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)) - a*c**2*x*sqrt(c*x)/(2*b*(a + b*x**2)**(sympy.S(1)/4)) + c*(c*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/4)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_939():
    f = sqrt(c*x)/(a + b*x**2)**(sympy.S(1)/4)
    F = sqrt(a)*sqrt(c*x)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(b)*(a + b*x**2)**(sympy.S(1)/4)) + x*sqrt(c*x)/(a + b*x**2)**(sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_940():
    f = 1/((c*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4))
    F = -2/(c*sqrt(c*x)*(a + b*x**2)**(sympy.S(1)/4)) + 2*sqrt(b)*sqrt(c*x)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(a)*c**2*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_941():
    f = 1/((c*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(1)/4))
    F = 4*b/(5*a*c**3*sqrt(c*x)*(a + b*x**2)**(sympy.S(1)/4)) - 2*(a + b*x**2)**(sympy.S(3)/4)/(5*a*c*(c*x)**(sympy.S(5)/2)) - 4*b**(sympy.S(3)/2)*sqrt(c*x)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(5*a**(sympy.S(3)/2)*c**4*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_942():
    f = 1/((c*x)**(sympy.S(11)/2)*(a + b*x**2)**(sympy.S(1)/4))
    F = -2*(a + b*x**2)**(sympy.S(3)/4)/(9*a*c*(c*x)**(sympy.S(9)/2)) - 8*b**2/(15*a**2*c**5*sqrt(c*x)*(a + b*x**2)**(sympy.S(1)/4)) + 4*b*(a + b*x**2)**(sympy.S(3)/4)/(15*a**2*c**3*(c*x)**(sympy.S(5)/2)) + 8*b**(sympy.S(5)/2)*sqrt(c*x)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(15*a**(sympy.S(5)/2)*c**6*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_943():
    f = (c*x)**(sympy.S(3)/2)/(a - b*x**2)**(sympy.S(1)/4)
    F = -sqrt(2)*a*c**(sympy.S(3)/2)*log(-sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(a - b*x**2)**(sympy.S(1)/4) + sqrt(b)*sqrt(c)*x/sqrt(a - b*x**2) + sqrt(c))/(16*b**(sympy.S(5)/4)) + sqrt(2)*a*c**(sympy.S(3)/2)*log(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(a - b*x**2)**(sympy.S(1)/4) + sqrt(b)*sqrt(c)*x/sqrt(a - b*x**2) + sqrt(c))/(16*b**(sympy.S(5)/4)) + sqrt(2)*a*c**(sympy.S(3)/2)*atan(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a - b*x**2)**(sympy.S(1)/4)) - 1)/(8*b**(sympy.S(5)/4)) + sqrt(2)*a*c**(sympy.S(3)/2)*atan(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a - b*x**2)**(sympy.S(1)/4)) + 1)/(8*b**(sympy.S(5)/4)) - c*sqrt(c*x)*(a - b*x**2)**(sympy.S(3)/4)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_944():
    f = 1/(sqrt(c*x)*(a - b*x**2)**(sympy.S(1)/4))
    F = -sqrt(2)*log(-sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(a - b*x**2)**(sympy.S(1)/4) + sqrt(b)*sqrt(c)*x/sqrt(a - b*x**2) + sqrt(c))/(4*b**(sympy.S(1)/4)*sqrt(c)) + sqrt(2)*log(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(a - b*x**2)**(sympy.S(1)/4) + sqrt(b)*sqrt(c)*x/sqrt(a - b*x**2) + sqrt(c))/(4*b**(sympy.S(1)/4)*sqrt(c)) + sqrt(2)*atan(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a - b*x**2)**(sympy.S(1)/4)) - 1)/(2*b**(sympy.S(1)/4)*sqrt(c)) + sqrt(2)*atan(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a - b*x**2)**(sympy.S(1)/4)) + 1)/(2*b**(sympy.S(1)/4)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_945():
    f = 1/((c*x)**(sympy.S(5)/2)*(a - b*x**2)**(sympy.S(1)/4))
    F = -2*(a - b*x**2)**(sympy.S(3)/4)/(3*a*c*(c*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_946():
    f = 1/((c*x)**(sympy.S(9)/2)*(a - b*x**2)**(sympy.S(1)/4))
    F = -2*(a - b*x**2)**(sympy.S(3)/4)/(3*a*c*(c*x)**(sympy.S(7)/2)) + 8*(a - b*x**2)**(sympy.S(7)/4)/(21*a**2*c*(c*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_947():
    f = 1/((c*x)**(sympy.S(13)/2)*(a - b*x**2)**(sympy.S(1)/4))
    F = -2*(a - b*x**2)**(sympy.S(3)/4)/(3*a*c*(c*x)**(sympy.S(11)/2)) + 16*(a - b*x**2)**(sympy.S(7)/4)/(21*a**2*c*(c*x)**(sympy.S(11)/2)) - 64*(a - b*x**2)**(sympy.S(11)/4)/(231*a**3*c*(c*x)**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_948():
    f = (c*x)**(sympy.S(5)/2)/(a - b*x**2)**(sympy.S(1)/4)
    F = a**(sympy.S(3)/2)*c**2*sqrt(c*x)*(-a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acsc(sqrt(b)*x/sqrt(a))/2, 2)/(2*b**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(1)/4)) - a*c**3*(a - b*x**2)**(sympy.S(3)/4)/(2*b**2*sqrt(c*x)) - c*(c*x)**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(3)/4)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_949():
    f = sqrt(c*x)/(a - b*x**2)**(sympy.S(1)/4)
    F = sqrt(a)*sqrt(c*x)*(-a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acsc(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(b)*(a - b*x**2)**(sympy.S(1)/4)) - c*(a - b*x**2)**(sympy.S(3)/4)/(b*sqrt(c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_950():
    f = 1/((c*x)**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(1)/4))
    F = -2*sqrt(b)*sqrt(c*x)*(-a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acsc(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(a)*c**2*(a - b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_951():
    f = 1/((c*x)**(sympy.S(7)/2)*(a - b*x**2)**(sympy.S(1)/4))
    F = -2*(a - b*x**2)**(sympy.S(3)/4)/(5*a*c*(c*x)**(sympy.S(5)/2)) - 4*b**(sympy.S(3)/2)*sqrt(c*x)*(-a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acsc(sqrt(b)*x/sqrt(a))/2, 2)/(5*a**(sympy.S(3)/2)*c**4*(a - b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_952():
    f = 1/((c*x)**(sympy.S(11)/2)*(a - b*x**2)**(sympy.S(1)/4))
    F = -2*(a - b*x**2)**(sympy.S(3)/4)/(9*a*c*(c*x)**(sympy.S(9)/2)) - 4*b*(a - b*x**2)**(sympy.S(3)/4)/(15*a**2*c**3*(c*x)**(sympy.S(5)/2)) - 8*b**(sympy.S(5)/2)*sqrt(c*x)*(-a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acsc(sqrt(b)*x/sqrt(a))/2, 2)/(15*a**(sympy.S(5)/2)*c**6*(a - b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_953():
    f = (c*x)**(sympy.S(3)/2)/(a + b*x**2)**(sympy.S(3)/4)
    F = sqrt(a)*(c*x)**(sympy.S(3)/2)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(b)*(a + b*x**2)**(sympy.S(3)/4)) + c*sqrt(c*x)*(a + b*x**2)**(sympy.S(1)/4)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_954():
    f = 1/(sqrt(c*x)*(a + b*x**2)**(sympy.S(3)/4))
    F = -2*sqrt(b)*(c*x)**(sympy.S(3)/2)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(a)*c**2*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_955():
    f = 1/((c*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/4))
    F = -2*(a + b*x**2)**(sympy.S(1)/4)/(3*a*c*(c*x)**(sympy.S(3)/2)) + 4*b**(sympy.S(3)/2)*(c*x)**(sympy.S(3)/2)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(3*a**(sympy.S(3)/2)*c**4*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_956():
    f = 1/((c*x)**(sympy.S(9)/2)*(a + b*x**2)**(sympy.S(3)/4))
    F = -2*(a + b*x**2)**(sympy.S(1)/4)/(7*a*c*(c*x)**(sympy.S(7)/2)) + 4*b*(a + b*x**2)**(sympy.S(1)/4)/(7*a**2*c**3*(c*x)**(sympy.S(3)/2)) - 8*b**(sympy.S(5)/2)*(c*x)**(sympy.S(3)/2)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(7*a**(sympy.S(5)/2)*c**6*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_957():
    f = 1/((c*x)**(sympy.S(13)/2)*(a + b*x**2)**(sympy.S(3)/4))
    F = -2*(a + b*x**2)**(sympy.S(1)/4)/(11*a*c*(c*x)**(sympy.S(11)/2)) + 20*b*(a + b*x**2)**(sympy.S(1)/4)/(77*a**2*c**3*(c*x)**(sympy.S(7)/2)) - 40*b**2*(a + b*x**2)**(sympy.S(1)/4)/(77*a**3*c**5*(c*x)**(sympy.S(3)/2)) + 80*b**(sympy.S(7)/2)*(c*x)**(sympy.S(3)/2)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(77*a**(sympy.S(7)/2)*c**8*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_958():
    f = (c*x)**(sympy.S(5)/2)/(a + b*x**2)**(sympy.S(3)/4)
    F = 3*a*c**(sympy.S(5)/2)*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(7)/4)) - 3*a*c**(sympy.S(5)/2)*atanh(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(7)/4)) + c*(c*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_959():
    f = sqrt(c*x)/(a + b*x**2)**(sympy.S(3)/4)
    F = -sqrt(c)*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/b**(sympy.S(3)/4) + sqrt(c)*atanh(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/b**(sympy.S(3)/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_960():
    f = 1/((c*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/4))
    F = -2*(a + b*x**2)**(sympy.S(1)/4)/(a*c*sqrt(c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_961():
    f = 1/((c*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(3)/4))
    F = -2*(a + b*x**2)**(sympy.S(1)/4)/(a*c*(c*x)**(sympy.S(5)/2)) + 8*(a + b*x**2)**(sympy.S(5)/4)/(5*a**2*c*(c*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_962():
    f = 1/((c*x)**(sympy.S(11)/2)*(a + b*x**2)**(sympy.S(3)/4))
    F = -2*(a + b*x**2)**(sympy.S(1)/4)/(a*c*(c*x)**(sympy.S(9)/2)) + 16*(a + b*x**2)**(sympy.S(5)/4)/(5*a**2*c*(c*x)**(sympy.S(9)/2)) - 64*(a + b*x**2)**(sympy.S(9)/4)/(45*a**3*c*(c*x)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_963():
    f = (c*x)**(sympy.S(3)/2)/(a - b*x**2)**(sympy.S(3)/4)
    F = -sqrt(a)*(c*x)**(sympy.S(3)/2)*(-a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acsc(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(b)*(a - b*x**2)**(sympy.S(3)/4)) - c*sqrt(c*x)*(a - b*x**2)**(sympy.S(1)/4)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_964():
    f = 1/(sqrt(c*x)*(a - b*x**2)**(sympy.S(3)/4))
    F = -2*sqrt(b)*(c*x)**(sympy.S(3)/2)*(-a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acsc(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(a)*c**2*(a - b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_965():
    f = 1/((c*x)**(sympy.S(5)/2)*(a - b*x**2)**(sympy.S(3)/4))
    F = -2*(a - b*x**2)**(sympy.S(1)/4)/(3*a*c*(c*x)**(sympy.S(3)/2)) - 4*b**(sympy.S(3)/2)*(c*x)**(sympy.S(3)/2)*(-a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acsc(sqrt(b)*x/sqrt(a))/2, 2)/(3*a**(sympy.S(3)/2)*c**4*(a - b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_966():
    f = 1/((c*x)**(sympy.S(9)/2)*(a - b*x**2)**(sympy.S(3)/4))
    F = -2*(a - b*x**2)**(sympy.S(1)/4)/(7*a*c*(c*x)**(sympy.S(7)/2)) - 4*b*(a - b*x**2)**(sympy.S(1)/4)/(7*a**2*c**3*(c*x)**(sympy.S(3)/2)) - 8*b**(sympy.S(5)/2)*(c*x)**(sympy.S(3)/2)*(-a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acsc(sqrt(b)*x/sqrt(a))/2, 2)/(7*a**(sympy.S(5)/2)*c**6*(a - b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_967():
    f = 1/((c*x)**(sympy.S(13)/2)*(a - b*x**2)**(sympy.S(3)/4))
    F = -2*(a - b*x**2)**(sympy.S(1)/4)/(11*a*c*(c*x)**(sympy.S(11)/2)) - 20*b*(a - b*x**2)**(sympy.S(1)/4)/(77*a**2*c**3*(c*x)**(sympy.S(7)/2)) - 40*b**2*(a - b*x**2)**(sympy.S(1)/4)/(77*a**3*c**5*(c*x)**(sympy.S(3)/2)) - 80*b**(sympy.S(7)/2)*(c*x)**(sympy.S(3)/2)*(-a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acsc(sqrt(b)*x/sqrt(a))/2, 2)/(77*a**(sympy.S(7)/2)*c**8*(a - b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_968():
    f = (c*x)**(sympy.S(5)/2)/(a - b*x**2)**(sympy.S(3)/4)
    F = 3*sqrt(2)*a*c**(sympy.S(5)/2)*log(-sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(a - b*x**2)**(sympy.S(1)/4) + sqrt(b)*sqrt(c)*x/sqrt(a - b*x**2) + sqrt(c))/(16*b**(sympy.S(7)/4)) - 3*sqrt(2)*a*c**(sympy.S(5)/2)*log(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(a - b*x**2)**(sympy.S(1)/4) + sqrt(b)*sqrt(c)*x/sqrt(a - b*x**2) + sqrt(c))/(16*b**(sympy.S(7)/4)) + 3*sqrt(2)*a*c**(sympy.S(5)/2)*atan(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a - b*x**2)**(sympy.S(1)/4)) - 1)/(8*b**(sympy.S(7)/4)) + 3*sqrt(2)*a*c**(sympy.S(5)/2)*atan(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a - b*x**2)**(sympy.S(1)/4)) + 1)/(8*b**(sympy.S(7)/4)) - c*(c*x)**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(1)/4)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_969():
    f = sqrt(c*x)/(a - b*x**2)**(sympy.S(3)/4)
    F = sqrt(2)*sqrt(c)*log(-sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(a - b*x**2)**(sympy.S(1)/4) + sqrt(b)*sqrt(c)*x/sqrt(a - b*x**2) + sqrt(c))/(4*b**(sympy.S(3)/4)) - sqrt(2)*sqrt(c)*log(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(a - b*x**2)**(sympy.S(1)/4) + sqrt(b)*sqrt(c)*x/sqrt(a - b*x**2) + sqrt(c))/(4*b**(sympy.S(3)/4)) + sqrt(2)*sqrt(c)*atan(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a - b*x**2)**(sympy.S(1)/4)) - 1)/(2*b**(sympy.S(3)/4)) + sqrt(2)*sqrt(c)*atan(sqrt(2)*b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a - b*x**2)**(sympy.S(1)/4)) + 1)/(2*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_970():
    f = 1/((c*x)**(sympy.S(3)/2)*(a - b*x**2)**(sympy.S(3)/4))
    F = -2*(a - b*x**2)**(sympy.S(1)/4)/(a*c*sqrt(c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_971():
    f = 1/((c*x)**(sympy.S(7)/2)*(a - b*x**2)**(sympy.S(3)/4))
    F = -2*(a - b*x**2)**(sympy.S(1)/4)/(a*c*(c*x)**(sympy.S(5)/2)) + 8*(a - b*x**2)**(sympy.S(5)/4)/(5*a**2*c*(c*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_972():
    f = 1/((c*x)**(sympy.S(11)/2)*(a - b*x**2)**(sympy.S(3)/4))
    F = -2*(a - b*x**2)**(sympy.S(1)/4)/(a*c*(c*x)**(sympy.S(9)/2)) + 16*(a - b*x**2)**(sympy.S(5)/4)/(5*a**2*c*(c*x)**(sympy.S(9)/2)) - 64*(a - b*x**2)**(sympy.S(9)/4)/(45*a**3*c*(c*x)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_973():
    f = (c*x)**(sympy.S(7)/2)/(a + b*x**2)**(sympy.S(5)/4)
    F = 5*a*c**3*sqrt(c*x)/(2*b**2*(a + b*x**2)**(sympy.S(1)/4)) - 5*a*c**(sympy.S(7)/2)*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(9)/4)) - 5*a*c**(sympy.S(7)/2)*atanh(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(9)/4)) + c*(c*x)**(sympy.S(5)/2)/(2*b*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_974():
    f = (c*x)**(sympy.S(3)/2)/(a + b*x**2)**(sympy.S(5)/4)
    F = -2*c*sqrt(c*x)/(b*(a + b*x**2)**(sympy.S(1)/4)) + c**(sympy.S(3)/2)*atan(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/b**(sympy.S(5)/4) + c**(sympy.S(3)/2)*atanh(b**(sympy.S(1)/4)*sqrt(c*x)/(sqrt(c)*(a + b*x**2)**(sympy.S(1)/4)))/b**(sympy.S(5)/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_975():
    f = 1/(sqrt(c*x)*(a + b*x**2)**(sympy.S(5)/4))
    F = 2*sqrt(c*x)/(a*c*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_976():
    f = 1/((c*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(5)/4))
    F = 2/(a*c*(c*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)) - 8*(a + b*x**2)**(sympy.S(3)/4)/(3*a**2*c*(c*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_977():
    f = 1/((c*x)**(sympy.S(9)/2)*(a + b*x**2)**(sympy.S(5)/4))
    F = 2/(a*c*(c*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(1)/4)) - 16*(a + b*x**2)**(sympy.S(3)/4)/(3*a**2*c*(c*x)**(sympy.S(7)/2)) + 64*(a + b*x**2)**(sympy.S(7)/4)/(21*a**3*c*(c*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_978():
    f = 1/((c*x)**(sympy.S(13)/2)*(a + b*x**2)**(sympy.S(5)/4))
    F = 2/(a*c*(c*x)**(sympy.S(11)/2)*(a + b*x**2)**(sympy.S(1)/4)) - 8*(a + b*x**2)**(sympy.S(3)/4)/(a**2*c*(c*x)**(sympy.S(11)/2)) + 64*(a + b*x**2)**(sympy.S(7)/4)/(7*a**3*c*(c*x)**(sympy.S(11)/2)) - 256*(a + b*x**2)**(sympy.S(11)/4)/(77*a**4*c*(c*x)**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_979():
    f = (c*x)**(sympy.S(13)/2)/(a + b*x**2)**(sympy.S(5)/4)
    F = 77*a**(sympy.S(5)/2)*c**6*sqrt(c*x)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(20*b**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(1)/4)) + 77*a**2*c**5*(c*x)**(sympy.S(3)/2)/(60*b**3*(a + b*x**2)**(sympy.S(1)/4)) - 11*a*c**3*(c*x)**(sympy.S(7)/2)/(30*b**2*(a + b*x**2)**(sympy.S(1)/4)) + c*(c*x)**(sympy.S(11)/2)/(5*b*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_980():
    f = (c*x)**(sympy.S(9)/2)/(a + b*x**2)**(sympy.S(5)/4)
    F = -7*a**(sympy.S(3)/2)*c**4*sqrt(c*x)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(2*b**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)) - 7*a*c**3*(c*x)**(sympy.S(3)/2)/(6*b**2*(a + b*x**2)**(sympy.S(1)/4)) + c*(c*x)**(sympy.S(7)/2)/(3*b*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_981():
    f = (c*x)**(sympy.S(5)/2)/(a + b*x**2)**(sympy.S(5)/4)
    F = 3*sqrt(a)*c**2*sqrt(c*x)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(b**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)) + c*(c*x)**(sympy.S(3)/2)/(b*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_982():
    f = sqrt(c*x)/(a + b*x**2)**(sympy.S(5)/4)
    F = -2*sqrt(c*x)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(a)*sqrt(b)*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_983():
    f = 1/((c*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(5)/4))
    F = -2/(a*c*sqrt(c*x)*(a + b*x**2)**(sympy.S(1)/4)) + 4*sqrt(b)*sqrt(c*x)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(a**(sympy.S(3)/2)*c**2*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_984():
    f = 1/((c*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(5)/4))
    F = -2/(5*a*c*(c*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)) + 12*b/(5*a**2*c**3*sqrt(c*x)*(a + b*x**2)**(sympy.S(1)/4)) - 24*b**(sympy.S(3)/2)*sqrt(c*x)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(5*a**(sympy.S(5)/2)*c**4*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_985():
    f = 1/((c*x)**(sympy.S(11)/2)*(a + b*x**2)**(sympy.S(5)/4))
    F = -2/(9*a*c*(c*x)**(sympy.S(9)/2)*(a + b*x**2)**(sympy.S(1)/4)) + 4*b/(9*a**2*c**3*(c*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)) - 8*b**2/(3*a**3*c**5*sqrt(c*x)*(a + b*x**2)**(sympy.S(1)/4)) + 16*b**(sympy.S(5)/2)*sqrt(c*x)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(3*a**(sympy.S(7)/2)*c**6*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_986():
    f = (c*x)**(sympy.S(5)/4)/(a + b*x**2)**(sympy.S(1)/4)
    F = 4*(c*x)**(sympy.S(9)/4)*(1 + b*x**2/a)**(sympy.S(1)/4)*hyper((sympy.S(1)/4, sympy.S(9)/8), (sympy.S(17)/8,), -b*x**2/a)/(9*c*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_987():
    f = (c*x)**(sympy.S(3)/4)/(a + b*x**2)**(sympy.S(1)/4)
    F = 4*(c*x)**(sympy.S(7)/4)*(1 + b*x**2/a)**(sympy.S(1)/4)*hyper((sympy.S(1)/4, sympy.S(7)/8), (sympy.S(15)/8,), -b*x**2/a)/(7*c*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_988():
    f = (c*x)**(sympy.S(1)/4)/(a + b*x**2)**(sympy.S(1)/4)
    F = 4*(c*x)**(sympy.S(5)/4)*(1 + b*x**2/a)**(sympy.S(1)/4)*hyper((sympy.S(1)/4, sympy.S(5)/8), (sympy.S(13)/8,), -b*x**2/a)/(5*c*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_989():
    f = 1/((c*x)**(sympy.S(1)/4)*(a + b*x**2)**(sympy.S(1)/4))
    F = 4*(c*x)**(sympy.S(3)/4)*(1 + b*x**2/a)**(sympy.S(1)/4)*hyper((sympy.S(1)/4, sympy.S(3)/8), (sympy.S(11)/8,), -b*x**2/a)/(3*c*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_990():
    f = 1/((c*x)**(sympy.S(3)/4)*(a + b*x**2)**(sympy.S(1)/4))
    F = 4*(c*x)**(sympy.S(1)/4)*(1 + b*x**2/a)**(sympy.S(1)/4)*hyper((sympy.S(1)/8, sympy.S(1)/4), (sympy.S(9)/8,), -b*x**2/a)/(c*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_991():
    f = 1/((c*x)**(sympy.S(5)/4)*(a + b*x**2)**(sympy.S(1)/4))
    F = -4*(1 + b*x**2/a)**(sympy.S(1)/4)*hyper((sympy.S(-1)/8, sympy.S(1)/4), (sympy.S(7)/8,), -b*x**2/a)/(c*(c*x)**(sympy.S(1)/4)*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_992():
    f = (c*x)**(sympy.S(5)/4)/(a + b*x**2)**(sympy.S(7)/4)
    F = 4*(c*x)**(sympy.S(9)/4)*(1 + b*x**2/a)**(sympy.S(3)/4)*hyper((sympy.S(9)/8, sympy.S(7)/4), (sympy.S(17)/8,), -b*x**2/a)/(9*a*c*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_993():
    f = (c*x)**(sympy.S(3)/4)/(a + b*x**2)**(sympy.S(7)/4)
    F = 4*(c*x)**(sympy.S(7)/4)*(1 + b*x**2/a)**(sympy.S(3)/4)*hyper((sympy.S(7)/8, sympy.S(7)/4), (sympy.S(15)/8,), -b*x**2/a)/(7*a*c*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_994():
    f = (c*x)**(sympy.S(1)/4)/(a + b*x**2)**(sympy.S(7)/4)
    F = 4*(c*x)**(sympy.S(5)/4)*(1 + b*x**2/a)**(sympy.S(3)/4)*hyper((sympy.S(5)/8, sympy.S(7)/4), (sympy.S(13)/8,), -b*x**2/a)/(5*a*c*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_995():
    f = 1/((c*x)**(sympy.S(1)/4)*(a + b*x**2)**(sympy.S(7)/4))
    F = 4*(c*x)**(sympy.S(3)/4)*(1 + b*x**2/a)**(sympy.S(3)/4)*hyper((sympy.S(3)/8, sympy.S(7)/4), (sympy.S(11)/8,), -b*x**2/a)/(3*a*c*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_996():
    f = 1/((c*x)**(sympy.S(3)/4)*(a + b*x**2)**(sympy.S(7)/4))
    F = 4*(c*x)**(sympy.S(1)/4)*(1 + b*x**2/a)**(sympy.S(3)/4)*hyper((sympy.S(1)/8, sympy.S(7)/4), (sympy.S(9)/8,), -b*x**2/a)/(a*c*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_997():
    f = 1/((c*x)**(sympy.S(5)/4)*(a + b*x**2)**(sympy.S(7)/4))
    F = -4*(1 + b*x**2/a)**(sympy.S(3)/4)*hyper((sympy.S(-1)/8, sympy.S(7)/4), (sympy.S(7)/8,), -b*x**2/a)/(a*c*(c*x)**(sympy.S(1)/4)*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_998():
    f = x**6*(a + b*x**2)**(sympy.S(1)/6)
    F = -81*3**(sympy.S(3)/4)*a**4*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*(a + b*x**2)**(sympy.S(1)/6)*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(2816*b**4*x*(a/(a + b*x**2))**(sympy.S(1)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)) + 81*a**3*x*(a + b*x**2)**(sympy.S(1)/6)/(2816*b**3) - 9*a**2*x**3*(a + b*x**2)**(sympy.S(1)/6)/(704*b**2) + 3*a*x**5*(a + b*x**2)**(sympy.S(1)/6)/(352*b) + 3*x**7*(a + b*x**2)**(sympy.S(1)/6)/22
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_999():
    f = x**4*(a + b*x**2)**(sympy.S(1)/6)
    F = 27*3**(sympy.S(3)/4)*a**3*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*(a + b*x**2)**(sympy.S(1)/6)*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(640*b**3*x*(a/(a + b*x**2))**(sympy.S(1)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)) - 27*a**2*x*(a + b*x**2)**(sympy.S(1)/6)/(640*b**2) + 3*a*x**3*(a + b*x**2)**(sympy.S(1)/6)/(160*b) + 3*x**5*(a + b*x**2)**(sympy.S(1)/6)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1000():
    f = x**2*(a + b*x**2)**(sympy.S(1)/6)
    F = -3*3**(sympy.S(3)/4)*a**2*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*(a + b*x**2)**(sympy.S(1)/6)*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(40*b**2*x*(a/(a + b*x**2))**(sympy.S(1)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)) + 3*a*x*(a + b*x**2)**(sympy.S(1)/6)/(40*b) + 3*x**3*(a + b*x**2)**(sympy.S(1)/6)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1001():
    f = (a + b*x**2)**(sympy.S(1)/6)
    F = 3**(sympy.S(3)/4)*a*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*(a + b*x**2)**(sympy.S(1)/6)*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(4*b*x*(a/(a + b*x**2))**(sympy.S(1)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)) + 3*x*(a + b*x**2)**(sympy.S(1)/6)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1002():
    f = (a + b*x**2)**(sympy.S(1)/6)/x**2
    F = -(a + b*x**2)**(sympy.S(1)/6)/x + 3**(sympy.S(3)/4)*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*(a + b*x**2)**(sympy.S(1)/6)*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*x*(a/(a + b*x**2))**(sympy.S(1)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1003():
    f = (a + b*x**2)**(sympy.S(1)/6)/x**4
    F = -(a + b*x**2)**(sympy.S(1)/6)/(3*x**3) - b*(a + b*x**2)**(sympy.S(1)/6)/(9*a*x) - 2*3**(sympy.S(3)/4)*b*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*(a + b*x**2)**(sympy.S(1)/6)*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(27*a*x*(a/(a + b*x**2))**(sympy.S(1)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1004():
    f = (a + b*x**2)**(sympy.S(1)/6)/x**6
    F = -(a + b*x**2)**(sympy.S(1)/6)/(5*x**5) - b*(a + b*x**2)**(sympy.S(1)/6)/(45*a*x**3) + 8*b**2*(a + b*x**2)**(sympy.S(1)/6)/(135*a**2*x) + 16*3**(sympy.S(3)/4)*b**2*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*(a + b*x**2)**(sympy.S(1)/6)*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(405*a**2*x*(a/(a + b*x**2))**(sympy.S(1)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1005():
    f = (a + b*x**2)**(sympy.S(1)/6)/x**8
    F = -(a + b*x**2)**(sympy.S(1)/6)/(7*x**7) - b*(a + b*x**2)**(sympy.S(1)/6)/(105*a*x**5) + 2*b**2*(a + b*x**2)**(sympy.S(1)/6)/(135*a**2*x**3) - 16*b**3*(a + b*x**2)**(sympy.S(1)/6)/(405*a**3*x) - 32*3**(sympy.S(3)/4)*b**3*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*(a + b*x**2)**(sympy.S(1)/6)*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(1215*a**3*x*(a/(a + b*x**2))**(sympy.S(1)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1006():
    f = x**6/(a + b*x**2)**(sympy.S(1)/6)
    F = -243*a**4*x/(896*b**3*(a/(a + b*x**2))**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(7)/6)*(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)) - 243*3**(sympy.S(1)/4)*a**4*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(1792*b**4*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) + 81*sqrt(2)*3**(sympy.S(3)/4)*a**4*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(896*b**4*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) - 243*a**3*x/(896*b**3*(a + b*x**2)**(sympy.S(1)/6)) + 81*a**2*x*(a + b*x**2)**(sympy.S(5)/6)/(448*b**3) - 9*a*x**3*(a + b*x**2)**(sympy.S(5)/6)/(56*b**2) + 3*x**5*(a + b*x**2)**(sympy.S(5)/6)/(20*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1007():
    f = x**4/(a + b*x**2)**(sympy.S(1)/6)
    F = 81*a**3*x/(224*b**2*(a/(a + b*x**2))**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(7)/6)*(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)) + 81*3**(sympy.S(1)/4)*a**3*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(448*b**3*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) - 27*sqrt(2)*3**(sympy.S(3)/4)*a**3*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(224*b**3*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) + 81*a**2*x/(224*b**2*(a + b*x**2)**(sympy.S(1)/6)) - 27*a*x*(a + b*x**2)**(sympy.S(5)/6)/(112*b**2) + 3*x**3*(a + b*x**2)**(sympy.S(5)/6)/(14*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1008():
    f = x**2/(a + b*x**2)**(sympy.S(1)/6)
    F = -9*a**2*x/(16*b*(a/(a + b*x**2))**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(7)/6)*(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)) - 9*3**(sympy.S(1)/4)*a**2*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(32*b**2*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) + 3*sqrt(2)*3**(sympy.S(3)/4)*a**2*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(16*b**2*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) - 9*a*x/(16*b*(a + b*x**2)**(sympy.S(1)/6)) + 3*x*(a + b*x**2)**(sympy.S(5)/6)/(8*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1009():
    f = (a + b*x**2)**(sympy.S(-1)/6)
    F = 3*a*x/(2*(a/(a + b*x**2))**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(7)/6)*(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)) + 3*3**(sympy.S(1)/4)*a*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(4*b*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) - sqrt(2)*3**(sympy.S(3)/4)*a*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(2*b*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) + 3*x/(2*(a + b*x**2)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1010():
    f = 1/(x**2*(a + b*x**2)**(sympy.S(1)/6))
    F = b*x/((a/(a + b*x**2))**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(7)/6)*(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)) + 3**(sympy.S(1)/4)*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(2*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) - sqrt(2)*3**(sympy.S(3)/4)*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) + b*x/(a*(a + b*x**2)**(sympy.S(1)/6)) - (a + b*x**2)**(sympy.S(5)/6)/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1011():
    f = 1/(x**4*(a + b*x**2)**(sympy.S(1)/6))
    F = -4*b**2*x/(9*a*(a/(a + b*x**2))**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(7)/6)*(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)) - 2*3**(sympy.S(1)/4)*b*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(9*a*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) + 4*sqrt(2)*3**(sympy.S(3)/4)*b*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(27*a*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) - (a + b*x**2)**(sympy.S(5)/6)/(3*a*x**3) - 4*b**2*x/(9*a**2*(a + b*x**2)**(sympy.S(1)/6)) + 4*b*(a + b*x**2)**(sympy.S(5)/6)/(9*a**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1012():
    f = 1/(x**6*(a + b*x**2)**(sympy.S(1)/6))
    F = -(a + b*x**2)**(sympy.S(5)/6)/(5*a*x**5) + 8*b**3*x/(27*a**2*(a/(a + b*x**2))**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(7)/6)*(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)) + 4*3**(sympy.S(1)/4)*b**2*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(27*a**2*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) - 8*sqrt(2)*3**(sympy.S(3)/4)*b**2*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(81*a**2*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) + 2*b*(a + b*x**2)**(sympy.S(5)/6)/(9*a**2*x**3) + 8*b**3*x/(27*a**3*(a + b*x**2)**(sympy.S(1)/6)) - 8*b**2*(a + b*x**2)**(sympy.S(5)/6)/(27*a**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1013():
    f = x**6/(a + b*x**2)**(sympy.S(5)/6)
    F = -81*3**(sympy.S(3)/4)*a**3*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*(a + b*x**2)**(sympy.S(1)/6)*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(128*b**4*x*(a/(a + b*x**2))**(sympy.S(1)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)) + 81*a**2*x*(a + b*x**2)**(sympy.S(1)/6)/(128*b**3) - 9*a*x**3*(a + b*x**2)**(sympy.S(1)/6)/(32*b**2) + 3*x**5*(a + b*x**2)**(sympy.S(1)/6)/(16*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1014():
    f = x**4/(a + b*x**2)**(sympy.S(5)/6)
    F = 27*3**(sympy.S(3)/4)*a**2*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*(a + b*x**2)**(sympy.S(1)/6)*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(40*b**3*x*(a/(a + b*x**2))**(sympy.S(1)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)) - 27*a*x*(a + b*x**2)**(sympy.S(1)/6)/(40*b**2) + 3*x**3*(a + b*x**2)**(sympy.S(1)/6)/(10*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1015():
    f = x**2/(a + b*x**2)**(sympy.S(5)/6)
    F = -3*3**(sympy.S(3)/4)*a*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*(a + b*x**2)**(sympy.S(1)/6)*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(4*b**2*x*(a/(a + b*x**2))**(sympy.S(1)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)) + 3*x*(a + b*x**2)**(sympy.S(1)/6)/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1016():
    f = (a + b*x**2)**(sympy.S(-5)/6)
    F = 3**(sympy.S(3)/4)*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*(a + b*x**2)**(sympy.S(1)/6)*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(b*x*(a/(a + b*x**2))**(sympy.S(1)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1017():
    f = 1/(x**2*(a + b*x**2)**(sympy.S(5)/6))
    F = -(a + b*x**2)**(sympy.S(1)/6)/(a*x) - 2*3**(sympy.S(3)/4)*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*(a + b*x**2)**(sympy.S(1)/6)*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*a*x*(a/(a + b*x**2))**(sympy.S(1)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1018():
    f = 1/(x**4*(a + b*x**2)**(sympy.S(5)/6))
    F = -(a + b*x**2)**(sympy.S(1)/6)/(3*a*x**3) + 8*b*(a + b*x**2)**(sympy.S(1)/6)/(9*a**2*x) + 16*3**(sympy.S(3)/4)*b*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*(a + b*x**2)**(sympy.S(1)/6)*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(27*a**2*x*(a/(a + b*x**2))**(sympy.S(1)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1019():
    f = 1/(x**6*(a + b*x**2)**(sympy.S(5)/6))
    F = -(a + b*x**2)**(sympy.S(1)/6)/(5*a*x**5) + 14*b*(a + b*x**2)**(sympy.S(1)/6)/(45*a**2*x**3) - 112*b**2*(a + b*x**2)**(sympy.S(1)/6)/(135*a**3*x) - 224*3**(sympy.S(3)/4)*b**2*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*(a + b*x**2)**(sympy.S(1)/6)*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(405*a**3*x*(a/(a + b*x**2))**(sympy.S(1)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1020():
    f = x**6/(a + b*x**2)**(sympy.S(7)/6)
    F = 1215*a**3*x/(224*b**3*(a/(a + b*x**2))**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(7)/6)*(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)) + 1215*3**(sympy.S(1)/4)*a**3*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(448*b**4*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) - 405*sqrt(2)*3**(sympy.S(3)/4)*a**3*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(224*b**4*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) + 1215*a**2*x/(224*b**3*(a + b*x**2)**(sympy.S(1)/6)) - 405*a*x*(a + b*x**2)**(sympy.S(5)/6)/(112*b**3) - 3*x**5/(b*(a + b*x**2)**(sympy.S(1)/6)) + 45*x**3*(a + b*x**2)**(sympy.S(5)/6)/(14*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1021():
    f = x**4/(a + b*x**2)**(sympy.S(7)/6)
    F = -81*a**2*x/(16*b**2*(a/(a + b*x**2))**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(7)/6)*(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)) - 81*3**(sympy.S(1)/4)*a**2*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(32*b**3*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) + 27*sqrt(2)*3**(sympy.S(3)/4)*a**2*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(16*b**3*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) - 81*a*x/(16*b**2*(a + b*x**2)**(sympy.S(1)/6)) - 3*x**3/(b*(a + b*x**2)**(sympy.S(1)/6)) + 27*x*(a + b*x**2)**(sympy.S(5)/6)/(8*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1022():
    f = x**2/(a + b*x**2)**(sympy.S(7)/6)
    F = 9*a*x/(2*b*(a/(a + b*x**2))**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(7)/6)*(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)) + 9*3**(sympy.S(1)/4)*a*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(4*b**2*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) - 3*sqrt(2)*3**(sympy.S(3)/4)*a*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(2*b**2*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) + 3*x/(2*b*(a + b*x**2)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1023():
    f = (a + b*x**2)**(sympy.S(-7)/6)
    F = -3*x/((a/(a + b*x**2))**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(7)/6)*(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)) - 3*3**(sympy.S(1)/4)*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(2*b*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(b*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1024():
    f = 1/(x**2*(a + b*x**2)**(sympy.S(7)/6))
    F = 4*b*x/(a*(a/(a + b*x**2))**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(7)/6)*(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)) + 3/(a*x*(a + b*x**2)**(sympy.S(1)/6)) + 2*3**(sympy.S(1)/4)*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(a*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) - 4*sqrt(2)*3**(sympy.S(3)/4)*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*a*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) + 4*b*x/(a**2*(a + b*x**2)**(sympy.S(1)/6)) - 4*(a + b*x**2)**(sympy.S(5)/6)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1025():
    f = 1/(x**4*(a + b*x**2)**(sympy.S(7)/6))
    F = 3/(a*x**3*(a + b*x**2)**(sympy.S(1)/6)) - 40*b**2*x/(9*a**2*(a/(a + b*x**2))**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(7)/6)*(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)) - 20*3**(sympy.S(1)/4)*b*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(9*a**2*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) + 40*sqrt(2)*3**(sympy.S(3)/4)*b*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(27*a**2*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) - 10*(a + b*x**2)**(sympy.S(5)/6)/(3*a**2*x**3) - 40*b**2*x/(9*a**3*(a + b*x**2)**(sympy.S(1)/6)) + 40*b*(a + b*x**2)**(sympy.S(5)/6)/(9*a**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1026():
    f = 1/(x**6*(a + b*x**2)**(sympy.S(7)/6))
    F = 3/(a*x**5*(a + b*x**2)**(sympy.S(1)/6)) - 16*(a + b*x**2)**(sympy.S(5)/6)/(5*a**2*x**5) + 128*b**3*x/(27*a**3*(a/(a + b*x**2))**(sympy.S(2)/3)*(a + b*x**2)**(sympy.S(7)/6)*(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)) + 64*3**(sympy.S(1)/4)*b**2*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(27*a**3*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) - 128*sqrt(2)*3**(sympy.S(3)/4)*b**2*sqrt(((a/(a + b*x**2))**(sympy.S(2)/3) + (a/(a + b*x**2))**(sympy.S(1)/3) + 1)/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (a/(a + b*x**2))**(sympy.S(1)/3))*elliptic_f(asin((-(a/(a + b*x**2))**(sympy.S(1)/3) + 1 + sqrt(3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(81*a**3*x*(a/(a + b*x**2))**(sympy.S(2)/3)*sqrt(-(1 - (a/(a + b*x**2))**(sympy.S(1)/3))/(-(a/(a + b*x**2))**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + b*x**2)**(sympy.S(1)/6)) + 32*b*(a + b*x**2)**(sympy.S(5)/6)/(9*a**3*x**3) + 128*b**3*x/(27*a**4*(a + b*x**2)**(sympy.S(1)/6)) - 128*b**2*(a + b*x**2)**(sympy.S(5)/6)/(27*a**4*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1027():
    f = x**7*(a + b*x**2)**p
    F = -a**3*(a + b*x**2)**(p + 1)/(2*b**4*(p + 1)) + 3*a**2*(a + b*x**2)**(p + 2)/(2*b**4*(p + 2)) - 3*a*(a + b*x**2)**(p + 3)/(2*b**4*(p + 3)) + (a + b*x**2)**(p + 4)/(2*b**4*(p + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1028():
    f = x**5*(a + b*x**2)**p
    F = a**2*(a + b*x**2)**(p + 1)/(2*b**3*(p + 1)) - a*(a + b*x**2)**(p + 2)/(b**3*(p + 2)) + (a + b*x**2)**(p + 3)/(2*b**3*(p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1029():
    f = x**3*(a + b*x**2)**p
    F = -a*(a + b*x**2)**(p + 1)/(2*b**2*(p + 1)) + (a + b*x**2)**(p + 2)/(2*b**2*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1030():
    f = x*(a + b*x**2)**p
    F = (a + b*x**2)**(p + 1)/(2*b*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1031():
    f = (a + b*x**2)**p/x
    F = -(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1032():
    f = (a + b*x**2)**p/x**3
    F = b*(a + b*x**2)**(p + 1)*hyper((2, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1033():
    f = (c*x)**m*(a + b*x**2)**p
    F = (c*x)**(m + 1)*(a + b*x**2)**p*hyper((-p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(c*(1 + b*x**2/a)**p*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1034():
    f = x**(-2*p - 7)*(a + b*x**2)**p
    F = -x**(-2*p - 6)*(a + b*x**2)**(p + 1)/(2*a*(p + 3)) + b*x**(-2*p - 4)*(a + b*x**2)**(p + 1)/(a**2*(p + 2)*(p + 3)) - b**2*x**(-2*p - 2)*(a + b*x**2)**(p + 1)/(a**3*(p + 1)*(p + 2)*(p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1035():
    f = x**(-2*p - 5)*(a + b*x**2)**p
    F = -x**(-2*p - 4)*(a + b*x**2)**(p + 1)/(2*a*(p + 2)) + b*x**(-2*p - 2)*(a + b*x**2)**(p + 1)/(2*a**2*(p + 1)*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_2_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1036():
    f = x**(-2*p - 3)*(a + b*x**2)**p
    F = -x**(-2*p - 2)*(a + b*x**2)**(p + 1)/(2*a*(p + 1))
    assert integrate(f, x) == F

