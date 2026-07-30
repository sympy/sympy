"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.6 Cosecant/4.6.4.2 (a+b csc)^m (d csc)^n (A+B csc+C csc^2).m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

A, B, C, a, b = symbols('A B C a b')

def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_4_2_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_plus_C_csc_pow_2_1():
    f = (a + b*csc(x))*(A + B*csc(x) + C*csc(x)**2)/sqrt(csc(x))
    F = -2*C*b*cos(x)*csc(x)**(sympy.S(3)/2)/3 + (-2*B*b - 2*C*a)*cos(x)*sqrt(csc(x)) - (2*B*b - 2*a*(A - C))*sqrt(sin(x))*sqrt(csc(x))*elliptic_e(x/2 - pi/4, 2) + (2*A*b + 2*B*a + 2*C*b/3)*sqrt(sin(x))*sqrt(csc(x))*elliptic_f(x/2 - pi/4, 2)
    assert integrate(f, x) == F

