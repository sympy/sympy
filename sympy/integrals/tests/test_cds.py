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
from sympy.integrals.cds import (cds_cancel_prim, cds_cancel_exp,
    cds_cancel_tan)
from sympy.utilities.pytest import raises

from sympy.abc import x, t, nu, z, a, y
t0, t1, t2 = symbols('t:3')
i = Symbol('i')

def test_cds_cancel_tan():
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1 + t**2, t)]})
    b0 = Poly(0, DE.t)
    b2 = 4*x
    c1 = -t**2 + 2*t - 8*x**2 + 1
    c2 = 2*(1 - 2*x)
    n = 2
    print cds_cancel_tan(b0, b2, c1, c2, DE, n)
