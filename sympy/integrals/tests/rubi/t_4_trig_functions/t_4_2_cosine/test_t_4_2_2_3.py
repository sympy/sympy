"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.2 Cosine/4.2.2.3 (g cos)^p (a+b cos)^m (c+d cos)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, c, e, f = symbols('a c e f')

def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_3_g_cos_pow_p_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_1():
    f = (a*cos(e + f*x) + a)**2*sec(e + f*x)**2/(c*cos(e + f*x) - c)
    F = -a**2*tan(e + f*x)/(c*f) - 3*a**2*atanh(sin(e + f*x))/(c*f) + 4*a**2*sin(e + f*x)/(c*f*(1 - cos(e + f*x)))
    assert integrate(f, x) == F

