"""Most of these tests come from the examples in Bronstein's book."""
from __future__ import with_statement
from sympy import (Poly, S, Function, log, symbols, exp, tan, sqrt,
    Symbol, Lambda, sin)
from sympy.integrals.risch import (gcdex_diophantine, frac_in, as_poly_1t,
    derivation, splitfactor, splitfactor_sqf, canonical_representation,
    hermite_reduce, polynomial_reduce, residue_reduce, residue_reduce_to_basic,
    integrate_primitive, integrate_hyperexponential_polynomial,
    integrate_hyperexponential, integrate_hypertangent_polynomial,
    integrate_nonlinear_no_specials, integer_powers, DifferentialExtension,
    risch_integrate, DecrementLevel, NonElementaryIntegral)
from sympy.integrals.cds import (cds_cancel_primitive, cds_cancel_exp,
    cds_cancel_tan, coupled_DE_system)
from sympy.utilities.pytest import raises

from sympy.abc import x, t, nu, z, a, y
t0, t1, t2 = symbols('t:3')
i = Symbol('i')

def test_cds_cancel_primitive():
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)],
        'L_K': [1], 'E_K': [], 'L_args': [x], 'E_args': []})
    b1 = Poly(2*x*t**3, DE.t)
    b2 = Poly(t + x, DE.t)
    c1 = Poly(t**2 + x**2, DE.t)
    c2 = Poly(2*t*x, DE.t)
    n = 4
    cds_cancel_primitive(Poly(sqrt(-1), DE.t), b1, b2, c1, c2, DE, n)

def test_cds_cancel_exp():
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)],
        'L_K': [1], 'E_K': [], 'L_args': [x], 'E_args': []})
    b1 = Poly(2*t**3, DE.t)
    b2 = Poly(4*t**2 + 2*t, DE.t)
    c1 = Poly(-t**2 + 2*t, DE.t)
    c2 = Poly(2*x*t, DE.t)
    n = 3
    cds_cancel_exp(Poly(sqrt(-1), DE.t), b1, b2, c1, c2, DE, n)
    b1 = Poly(t**3 + t**2 + 1, DE.t)
    b2 = Poly(t**4 + t**2 + 2*t + x**3, DE.t)
    c1 = Poly(-t**2 + 2*t + x**2, DE.t)
    c2 = Poly(2*x*t, DE.t)
    n = 2
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    cds_cancel_exp(Poly(sqrt(-1), DE.t), b1, b2, c1, c2, DE, n)
    b1 = Poly(t**2, DE.t)
    b2 = Poly(2*t**2, DE.t)
    c1 = Poly(-t**2 + 2*t , DE.t)
    c2 = Poly(2*t, DE.t)
    n = 2
    cds_cancel_exp(Poly(sqrt(-1), DE.t), b1, b2, c1, c2, DE, n)


def test_cds_cancel_tan():
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)]})
    b0 = Poly(5*t + 2*t**3, DE.t)
    b2 = Poly(4*t**2, DE.t)
    c1 = Poly(-t**2 + 2*t - 8*x**2 + 1, DE.t)
    c2 = Poly(2*(1 - 2*x), DE.t)
    n = 2
#    cds_cancel_tan(b0, b2, c1, c2, DE, n)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)],
       'L_K': [], 'E_K': [], 'L_args': [], 'E_args': []})
    b0 = Poly(0, DE.t)
    b2 = Poly(4*x, DE.t)
    c1 = Poly(8*x**2*1/(t**2 + 1) - t**2*1/(t**2 + 1) + 2*t*1/(t**2 + 1) + 1/(t**2 + 1), 1/(t**2 + 1))
    c2 = Poly(-4*x*1/(t**2 + 1) + 2*1/(t**2 + 1), 1/(t**2 + 1))
    n = 2
    cds_cancel_tan(b0, b2, c1, c2, DE, n)

def test_coupled_DE_system():
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)],
       'L_K': [], 'E_K': [1], 'L_args': [], 'E_args': [x]})
    g1 = Poly(2*t + 1, t)
    g2 = Poly(2*t**2 + 2*t + 1)
    assert coupled_DE_system(Poly(t + 1, t), Poly(t - 1, t), g1, g2, DE) == \
       (Poly(t, t), Poly(t + 1, t))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)],
       'L_K': [1], 'E_K': [], 'L_args': [x], 'E_args': []})
    g1 = Poly(t + 1 + x*t - x*t**3, t)
    g2 = Poly(2 + 2*t + 3*x*t**3 + 2*x*t)
    assert coupled_DE_system(Poly(t**2 + 1, t), Poly(t**2, t), g1, g2, DE) == \
       Poly((-x/(1 + I) + 3*I*x/(1 + I))*t, t)

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t1), Poly(1/(x + 1), t2)],
       'L_K': [1, 2], 'E_K': [], 'L_args': [x, x + 1], 'E_args': []})
    g1 = Poly(t + 1 +x*t - x*t**3, t)
    g2 = Poly(2 + 2*t + 3*x*t**3 + 2*x*t)
    assert coupled_DE_system(Poly(t**2 + 1, t), Poly(t**2, t), g1, g2, DE) == \
       Poly((-x/(1 + I) + 3*I*x/(1 + I))*t, t)


    g1 = Poly(t**3 - 3*t**2*x - 4*t**2 - t*x - t - x - 1 + 1/x, t)
    g2 = Poly(2*t**3 + 3*t**2 + t**2*x - t*x**2 - 2*t*x + 2*t + 2 + 2/x, t)
    assert coupled_DE_system(Poly(t**2 + 1, t), Poly(2*t + t*x, t), g1, g2, DE) == \
       (Poly(t - x, t), Poly(2*t + 1, t))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)],
       'L_K': [], 'E_K': [], 'L_args': [], 'E_args': []})
    b1 = Poly(5*t + 2*t**3 + x**2, DE.t)
    b2 = Poly(4*t**2 + x, DE.t)
    c1 = Poly(-t**2 + 2*t - 8*x**2 + 1, DE.t)
    c2 = Poly(2*(1 - 2*x), DE.t)
    n = 4
    coupled_DE_system(b1, b2, c1, c2, DE)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t1), Poly(1/(x + 1), t2)],
       'L_K': [1, 2], 'E_K': [], 'L_args': [x, x + 1], 'E_args': []})
    coupled_DE_system(Poly(t, t), Poly(1, t),Poly(t**2 - 3*t + 1, t), Poly(-t**2 + t + 1, t), DE)
