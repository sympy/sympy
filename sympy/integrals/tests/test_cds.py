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
    risch_integrate, DecrementLevel, NonElementaryIntegralException)
from sympy.integrals.cds import (cds_cancel_primitive, cds_cancel_exp,
    cds_cancel_tan, coupled_DE_system)
from sympy.utilities.pytest import raises

from sympy.abc import x, t, nu, z, a, y
t0, t1, t2 = symbols('t:3')
i = Symbol('i')

def test_cds_cancel_primitive():

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t1),
        Poly(2*x/(x**2 + 1), t2)], 'L_K': [1, 2], 'E_K': [], 'L_args': [x, x**2 + 1],
        'E_args': []})
    y1 = Poly(t2 + t1, t2)
    y2 = Poly(t2**2, t2)
    b1 = Poly(x, t2)
    b2 = Poly(2*x, t2)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, DE.t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, DE.t)
    n = 5
    assert cds_cancel_primitive(Poly(sqrt(-1), t2), b1, b2, a, b, DE, n) == \
        (y1, y2)
    assert cds_cancel_primitive(Poly(sqrt(-1), t), Poly(x, t), Poly(2*x, t),
        Poly(2*x/(x**2 + 1) + x*t2 - 2*x*t1), Poly(1/x + 2*x*t2 + x*t1), DE, 5) == \
        (t2, t1)
    y1 = Poly(t2 + t1, t2)
    y2 = Poly(t2**2, t2)
    b1 = Poly(x + 3*x**2, t2)
    b2 = Poly(2*x + 2*x, t2)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, DE.t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, DE.t)
    n = 3
    assert cds_cancel_primitive(Poly(sqrt(-1), t2), b1, b2, a, b, DE, n) == \
        (y1, y2)
    y1 = Poly(t2*t1**3, t2)
    y2 = Poly(t2**2 + t1, t2)
    b1 = Poly(x + 3*x**2, t2)
    b2 = Poly(x - 1, t2)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, DE.t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, DE.t)
    n = 4
    assert cds_cancel_primitive(Poly(sqrt(-1), t2), b1, b2, a, b, DE, n) == \
        (y1, y2)
    y1 = Poly(t2, t2)
    y2 = Poly(t2, t2)
    b1 = Poly(x + 3*x**2, t2)
    b2 = Poly(2*x + 2*x, t2)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, DE.t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, DE.t)
    n = 2
    assert cds_cancel_primitive(Poly(sqrt(-1), t2), b1, b2, a, b, DE, n) == \
        (y1, y2)
    y1 = Poly(t2**10, t2)
    y2 = Poly(t2, t2)
    b1 = Poly(x + 3*x**2, t2)
    b2 = Poly(2*x + 2*x, t2)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, DE.t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, DE.t)
    n = 11
    assert cds_cancel_primitive(Poly(sqrt(-1), t2), b1, b2, a, b, DE, n) == \
        (y1, y2)
    y1 = Poly(t2/x, t2)
    y2 = Poly(t2 + t1**10, t2)
    b1 = Poly(x + 3*x**2, t2)
    b2 = Poly(2*x + 2*x, t2)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, DE.t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, DE.t)
    n = 11
    assert cds_cancel_primitive(Poly(sqrt(-1), t2), b1, b2, a, b, DE, n) == \
        (y1, y2)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)],
        'L_K': [1], 'E_K': [], 'L_args': [x], 'E_args': [], 'Tfuncs': [log]})
    y1 = Poly(t + x**2, t)
    y2 = Poly(t, t)
    b1 = Poly(x, t)
    b2 = Poly(2*x, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, DE.t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, DE.t)
    n = 2
    assert cds_cancel_primitive(Poly(sqrt(-1), t2), b1, b2, a, b, DE, n) == \
        (y1, y2)
    y1 = Poly(t**2, t)
    y2 = Poly(t/x, t)
    b1 = Poly(x + 1, t)
    b2 = Poly(2*x, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, DE.t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, DE.t)
    n = 2
    assert cds_cancel_primitive(Poly(sqrt(-1), t), b1, b2, a, b, DE, n) == \
        (y1, y2)
    y1 = Poly(1 + t, t)
    y2 = Poly(t/x, t)
    b1 = Poly(x + 1, t)
    b2 = Poly(2*x, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, DE.t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, DE.t)
    n = 1
    assert cds_cancel_primitive(Poly(sqrt(-1), t), b1, b2, a, b, DE, n) == \
        (y1, y2)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(2/x, t)],
       'L_K': [1], 'E_K': [], 'L_args': [x], 'E_args': []})
    y1 = Poly(1 + t, t)
    y2 = Poly(t, t)
    b1 = Poly(x + 1, t)
    b2 = Poly(2*x, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, DE.t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, DE.t)
    n = 3
    assert cds_cancel_primitive(Poly(sqrt(-1), t), b1, b2, a, b, DE, n) == \
        (y1, y2)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t1, t1), Poly(t1/(t1 + 1), t2)],
       'L_K': [1], 'E_K': [], 'L_args': [x], 'E_args': []})
    y1 = Poly(1 + t2, t2)
    y2 = Poly(t2*t1, t2)
    b1 = Poly(x + 1, t2)
    b2 = Poly(2*x, t2)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, DE.t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, DE.t)
    n = 3
    assert cds_cancel_primitive(Poly(sqrt(-1), t), b1, b2, a, b, DE, n) == \
        (y1, y2)


def test_cds_cancel_exp():
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)],
        'L_K': [], 'E_K': [1], 'L_args': [], 'E_args': [x], 'Tfuncs': [exp]})
    y1 = Poly(t, t)
    y2 = Poly(t, t)
    b1 = Poly(x, t)
    b2 = Poly(2*x, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, DE.t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, DE.t)
    n = 7
    assert cds_cancel_exp(Poly(sqrt(-1), t), b1, b2, a, b, DE, n) == \
        (y1, y2)
    y1 = Poly(t**2, t)
    y2 = Poly(x, t)
    b1 = Poly(x/10, t)
    b2 = Poly(x**2, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, DE.t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, DE.t)
    n = 3
    assert cds_cancel_exp(Poly(sqrt(-1), t), b1, b2, a, b, DE, n) == \
        (y1, y2)
    y1 = Poly(t + x**2, t)
    y2 = Poly(t*x, t)
    b1 = Poly(x**2, t)
    b2 = Poly(2*x + 3, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, DE.t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, DE.t)
    n = 3
    assert cds_cancel_exp(Poly(sqrt(-1), t), b1, b2, a, b, DE, n) == \
        (y1, y2)
    y1 = Poly(t + 1, t)
    y2 = Poly(t**2, t)
    b1 = Poly(x + 4, t)
    b2 = Poly(2*x, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, DE.t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, DE.t)
    n = 4
    assert cds_cancel_exp(Poly(sqrt(-1), t), b1, b2, a, b, DE, n) == \
        (y1, y2)

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t1, t1), Poly(2*x*t2, t2)],
        'L_K': [], 'E_K': [1, 2], 'L_args': [], 'E_args': [x, x**2]})
    y1 = Poly(t2 + t1, t2)
    y2 = Poly(t2**2, t2)
    b1 = Poly(x + 4, t2)
    b2 = Poly(2*x**2 + 10, t2)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, t2)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, t2)
    n = 4
    assert cds_cancel_exp(Poly(sqrt(-1), t), b1, b2, a, b, DE, n) == \
        (y1, y2)
    y1 = Poly(t2*t1, t2)
    y2 = Poly(t2 + 3*t1, t2)
    b1 = Poly(x + x**2, t2)
    b2 = Poly(2*x**2 + 2*x + 1, t2)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, t2)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, t2)
    n = 4
    assert cds_cancel_exp(Poly(sqrt(-1), t), b1, b2, a, b, DE, n) == \
        (y1, y2)
    y1 = Poly(t2 + t1**2, t2)
    y2 = Poly(t2*t1, t2)
    b1 = Poly(x + x**2, t2)
    b2 = Poly(2*x**2 + 1, t2)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, t2)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, t2)
    n = 4
    assert cds_cancel_exp(Poly(sqrt(-1), t), b1, b2, a, b, DE, n) == \
        (y1, y2)

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(-t/x**2, t)],
        'L_K': [], 'E_K': [1], 'L_args': [], 'E_args': [1/x]})
    y1 = Poly(t + x**2, t)
    y2 = Poly(t**2 + x, t)
    b1 = Poly(x, t2)
    b2 = Poly(2*x + 1, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, t)
    n = 4
    assert cds_cancel_exp(Poly(sqrt(-1), t), b1, b2, a, b, DE, n) == \
        (y1, y2)
    y1 = Poly(t/x, t)
    y2 = Poly(t**2/x**2, t)
    b1 = Poly(x + 10, t2)
    b2 = Poly(2*x + x**2, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, t)
    n = 4
    assert cds_cancel_exp(Poly(sqrt(-1), t), b1, b2, a, b, DE, n) == \
        (y1, y2)
    y1 = Poly(t, t)
    y2 = Poly(t + x, t)
    b1 = Poly(x + x**3, t2)
    b2 = Poly(2*x + 1, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, t)
    n = 4
    assert cds_cancel_exp(Poly(sqrt(-1), t), b1, b2, a, b, DE, n) == \
        (y1, y2)


def test_cds_cancel_tan():
    #structure theorem required
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)],
        'L_K': [], 'E_K': [], 'T_K': [1], 'AT_K': [], 'L_args': [], 'E_args': [],
        'T_args': [x], 'AT_args': []})
    b0 = Poly(5*t + 2*t**3, DE.t)
    b2 = Poly(4*t**2, DE.t)
    c1 = Poly(-t**2 + 2*t - 8*x**2 + 1, DE.t)
    c2 = Poly(2*(1 - 2*x), DE.t)
    n = 2
    raises(NotImplementedError, lambda: cds_cancel_tan(b0, b2, c1, c2, DE, n))


def test_coupled_DE_system():
    DE = DifferentialExtension(extension={'D': [Poly(1, x)], 'L_K': [], 'E_K': [],
        'L_args': [],  'E_args': []})
    assert coupled_DE_system(Poly(x, x), Poly(2*x, x), Poly(-1/x**2),
        Poly(-1/(2*x**2) + S(5)/2), DE) == (1/x, 1/(2*x))
    raises(NonElementaryIntegralException, lambda: coupled_DE_system(Poly(x, x), Poly(2*x, x), Poly(-1/x**2),
        Poly(-1/(2*x**2)), DE))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)], 'L_K': [], 'E_K': [1],
        'L_args': [],  'E_args': [x]})
    y1 = Poly(t, t)
    y2 = Poly(t + x, t)
    b1 = Poly(x + x**3, t)
    b2 = Poly(2*x + 1, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, t)
    assert coupled_DE_system(b1, b2, Poly((x**3 - x)*t - 2*x**2 - x, t),
        Poly((x**3 + 3*x + 2)*t + x**4 + x**2 + 1, t), DE) == (y1, y2)
    y1 = Poly(t, t)
    y2 = Poly(t, t)
    b1 = Poly(x, t)
    b2 = Poly(2*x, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, t)
    assert coupled_DE_system(b1, b2, a, b, DE) == (y1, y2)
    y1 = Poly(t + t**2, t)
    y2 = Poly(t, t)
    b1 = Poly(3*x/2, t)
    b2 = Poly(x, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, t)
    assert coupled_DE_system(b1, b2, a, b, DE) == (y1, y2)

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)], 'L_K': [1], 'E_K': [],
        'L_args': [x],  'E_args': []})
    y1 = Poly(t + t**2, t)
    y2 = Poly(2*t + x, t)
    b1 = Poly(x + 5*x**2, t)
    b2 = Poly(2*x, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, t)
    assert coupled_DE_system(b1, b2, a, b, DE) == (y1, y2)
    y1 = Poly(t + x, t)
    y2 = Poly(t, t)
    b1 = Poly(3*x/2, t)
    b2 = Poly(x, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, t)
    assert coupled_DE_system(b1, b2, a, b, DE) == (y1, y2)
    y1 = Poly(t/x, t)
    y2 = Poly(x, t)
    b1 = Poly(x + 1, t)
    b2 = Poly(x/2, t)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, t)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, t)
    assert coupled_DE_system(b1, b2, a, b, DE) == (y1, y2)

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t1), Poly(t2, t2)], 'L_K': [1], 'E_K': [1],
        'L_args': [x],  'E_args': [x]})
    y1 = Poly(t1 + t2**2, t2)
    y2 = Poly(2*t1 + x, t2)
    b1 = Poly(x + 5*x**2, t2)
    b2 = Poly(2*x, t2)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, t2)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, t2)
    assert coupled_DE_system(b1, b2, a, b, DE) == (y1, y2)
    y1 = Poly(t1 + x, t2)
    y2 = Poly(t1, t2)
    b1 = Poly(3*x/2, t2)
    b2 = Poly(x, t2)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, t2)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, t2)
    assert coupled_DE_system(b1, b2, a, b, DE) == (y1, y2)
    y1 = Poly(t1/x, t2)
    y2 = Poly(x*t2, t2)
    b1 = Poly(x + 1, t2)
    b2 = Poly(x/2, t2)
    a = Poly(derivation(y1, DE) + b1*y1 - b2*y2, t2)
    b = Poly(derivation(y2, DE) + b2*y1 + b1*y2, t2)
    assert coupled_DE_system(b1, b2, a, b, DE) == (y1, y2)
