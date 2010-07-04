"""Most of these tests come from the examples in Bronstein's book."""
from sympy import Poly, S, Function, log
from sympy.integrals.risch import (gcdex_diophantine, derivation, splitfactor,
    splitfactor_sqf, canonical_representation, hermite_reduce,
    polynomial_reduce, residue_reduce, integrate_hypertangent_polynomial,
    integrate_nonlinear_no_specials,)
from sympy.utilities.pytest import XFAIL, skip

from sympy.abc import x, t, nu, z, a

def test_gcdex_diophantine():
    assert gcdex_diophantine(Poly(x**4 - 2*x**3 - 6*x**2 + 12*x + 15),
    Poly(x**3 + x**2 - 4*x - 4), Poly(x**2 - 1)) == \
        (Poly((-x**2 + 4*x - 3)/5), Poly((x**3 - 7*x**2 + 16*x - 10)/5))

def test_derivation():
    p = Poly(4*x**4*t**5 + (-4*x**3 - 4*x**4)*t**4 + (-3*x**2 + 2*x**3)*t**3 +
        (2*x + 7*x**2 + 2*x**3)*t**2 + (1 - 4*x - 4*x**2)*t - 1 + 2*x, t, domain='ZZ(x)')
    D = Poly(-t**2 - 3/(2*x)*t + 1/(2*x), t, domain='ZZ(x)')
    assert derivation(p, D, x, t) == Poly(-20*x**4*t**6 + (2*x**3 + 16*x**4)*t**5 +
        (21*x**2 + 12*x**3)*t**4 + (7*x/2 - 25*x**2 - 12*x**3)*t**3 +
        (-5 - 15*x/2 + 7*x**2)*t**2 - (3 - 8*x - 10*x**2 - 4*x**3)/(2*x)*t +
        (1 - 4*x**2)/(2*x), t, domain='ZZ(x)')
    assert derivation(Poly(1, t), D, x, t) == Poly(0, t)
    assert derivation(Poly(t, t), D, x, t) == D
    assert derivation(Poly(t**2 + 1/x*t + (1 - 2*x)/(4*x**2), t, domain='ZZ(x)'), D, x, t) == \
        Poly(-2*t**3 - 4/x*t**2 - (5 - 2*x)/(2*x**2)*t - (1 - 2*x)/(2*x**3), t, domain='ZZ(x)')

def test_splitfactor():
    p = Poly(4*x**4*t**5 + (-4*x**3 - 4*x**4)*t**4 + (-3*x**2 + 2*x**3)*t**3 +
        (2*x + 7*x**2 + 2*x**3)*t**2 + (1 - 4*x - 4*x**2)*t - 1 + 2*x, t, domain='ZZ(x)')
    D = Poly(-t**2 - 3/(2*x)*t + 1/(2*x), t, domain='ZZ(x)')
    assert splitfactor(p, D, x, t) == (Poly(4*x**4*t**3 + (-8*x**3 - 4*x**4)*t**2 +
        (4*x**2 + 8*x**3)*t - 4*x**2, t, domain='ZZ(x)'),
        Poly(t**2 + 1/x*t + (1 - 2*x)/(4*x**2), t, domain='ZZ(x)'))
    assert splitfactor(Poly(x, t), D, x, t) == (Poly(x, t), Poly(1, t))
    r = Poly(-4*x**4*z**2 + 4*x**6*z**2 - z*x**3 - 4*x**5*z**3 + 4*x**3*z**3 + x**4 + z*x**5 - x**6, t, domain='ZZ[x,z]')
    D = Poly(1/x, t, domain='ZZ(x)')
    assert splitfactor(r, D, x, t, coefficientD=True) == \
        (Poly(x*z - x**2 - z*x**3 + x**4, t, domain='ZZ[x,z]'),
        Poly(-x**2 + 4*x**2*z**2, t, domain='ZZ[x,z]'))
    assert splitfactor_sqf(r, D, x, t, coefficientD=True) == \
        (((Poly(x*z - x**2 - z*x**3 + x**4, t, domain='ZZ[x,z]'), 1),),
        ((Poly(-x**2 + 4*x**2*z**2, t, domain='ZZ[x,z]'), 1),))
    assert splitfactor(Poly(0, t), D, x, t) == (Poly(0, t), Poly(1, t))
    assert splitfactor_sqf(Poly(0, t), D, x, t) == (((Poly(0, t), 1),), ())

def test_canonical_representation():
    D = Poly(1 + t**2, t)
    assert canonical_representation(Poly(x - t, t), Poly(t**2, t), D, x, t) == \
        (Poly(0, t, domain='ZZ[x]'), (Poly(0, t, domain='ZZ[x]'),
        Poly(1, t, domain='ZZ')), (Poly(-t + x, t, domain='ZZ[x]'),
        Poly(t**2, t, domain='ZZ')))
    D = Poly(t**2 + 1, t)
    assert canonical_representation(Poly(t**5 + t**3 + x**2*t + 1, t),
    Poly((t**2 + 1)**3, t), Poly(1 + t**2, t), x, t) == \
        (Poly(0, t, domain='ZZ[x]'), (Poly(t**5 + t**3 + x**2*t + 1, t, domain='ZZ[x]'),
        Poly(t**6 + 3*t**4 + 3*t**2 + 1, t, domain='ZZ')), (Poly(0, t, domain='ZZ[x]'),
        Poly(1, t, domain='ZZ')))

def test_hermite_reduce():
    D = Poly(t**2 + 1, t)
    assert hermite_reduce(Poly(x - t, t), Poly(t**2, t), D, x, t) == \
        ((Poly(-x, t, domain='ZZ[x]'), Poly(t, t, domain='ZZ[x]')),
        (Poly(0, t, domain='ZZ[x]'), Poly(t, t, domain='ZZ')),
        (Poly(-x, t, domain='ZZ[x]'), Poly(1, t, domain='ZZ[x]')))
    D = Poly(-t**2 - t/x - (1 - nu**2/x**2), t)
    assert hermite_reduce(Poly(x**2*t**5 + x*t**4 - nu**2*t**3 - x*(x**2 + 1)*t**2 -
    (x**2 - nu**2)*t - x**5/4, t), Poly(x**2*t**4 + x**2*(x**2 + 2)*t**2 + x**2 +
    x**4 + x**6/4, t), D, x, t) == \
        ((Poly(-1 - x**2/4, t), Poly(t**2 + 1 + x**2/2, t)),
        (Poly((-2*nu**2 - x**4)/(2*x**2)*t - (1 + x**2)/x, t),
        Poly(t**2 + 1 + x**2/2, t)), (Poly(t + 1/x, t), Poly(1, t)))
    D = Poly(1/x, t)
    assert hermite_reduce(Poly(-t**2 + 2*t + 2, t),
    Poly(-x*t**2 + 2*x*t - x, t), D, x, t) == \
        ((Poly(3, t, domain='ZZ(x)'), Poly(t - 1, t, domain='ZZ(x)')),
        (Poly(0, t, domain='ZZ(x)'), Poly(t - 1, t, domain='ZZ(x)')),
        (Poly(1/x, t, domain='ZZ(x)'), Poly(1, t, domain='ZZ(x)')))
    assert hermite_reduce(Poly(-x**2*t**6 + (-1 - 2*x**3 + x**4)*t**3 +
    (-3 - 3*x**4)*t**2 - 2*x*t - x - 3*x**2, t, domain='ZZ[x]'),
    Poly(x**4*t**6 - 2*x**2*t**3 + 1, t, domain='ZZ[x]'), D, x, t) == \
        ((Poly(t + (1 + x**4)/x**3, t, domain='ZZ(x)'), Poly(t**3 - 1/x**2, t,
        domain='ZZ(x)')), (Poly(0, t, domain='ZZ(x)'), Poly(t**3 - 1/x**2, t,
        domain='ZZ(x)')), (Poly(-1/x**2, t, domain='ZZ(x)'), Poly(1, t, domain='ZZ(x)')))


def test_polynomial_reduce():
    D = Poly(1 + t**2, t)
    assert polynomial_reduce(Poly(1 + x*t + t**2, t), D, x, t) == \
        (Poly(t, t), Poly(x*t, t))
    assert polynomial_reduce(Poly(0, t), D, x, t) == \
        (Poly(0, t, domain='ZZ'), Poly(0, t, domain='ZZ'))

def test_residue_reduce():
    a = Poly(2*t**2 - t - x**2, t)
    d = Poly(t**3 - x**2*t, t)
    D = Poly(1/x, t)
    assert residue_reduce(a, d, D, x, t, z, invert=False) == \
        ([(Poly(z**2 - S(1)/4, z, domain='ZZ(x)'), Poly((1 + 3*x*z - 6*z**2 -
        2*x**2 + 4*x**2*z**2)*t - x*z + x**2 + 2*x**2*z**2 - 2*z*x**3, t,
        domain='ZZ[x,z]'))], False)
    assert residue_reduce(a, d, D, x, t, z, invert=True) == \
        ([(Poly(z**2 - S(1)/4, z, domain='ZZ(x)'), Poly(t + 2*x*z, t,
        domain='ZZ[x,z]'))], False)
    assert residue_reduce(Poly(-2/x, t, domain='ZZ(x)'), Poly(t**2 - 1, t,
    domain='ZZ(x)'), D, x, t, z, invert=False) == \
        ([(Poly(z**2 - 1, z, domain='ZZ(x)'), Poly(-z*t - 1, t,
        domain='ZZ(x,z)'))], True)
    assert residue_reduce(Poly(-2/x, t, domain='ZZ(x)'), Poly(t**2 - 1, t,
    domain='ZZ(x)'), D, x, t, z, invert=True) == \
        ([(Poly(z**2 - 1, z, domain='ZZ(x)'), Poly(t + z, t,
        domain='ZZ[z]'))], True)
    D = Poly(-t**2 - t/x - (1 - nu**2/x**2), t)
    assert residue_reduce(Poly((-2*nu**2 - x**4)/(2*x**2)*t - (1 + x**2)/x, t),
    Poly(t**2 + 1 + x**2/2, t), D, x, t, z) == \
        ([(Poly(z + S(1)/2, z, domain='QQ'), Poly(t**2 + 1 + x**2/2, t,
        domain='QQ[x]'))], True)
    D = Poly(1 + t**2, t)
    assert residue_reduce(Poly(-2*x*t + 1 - x**2, t, domain='ZZ(x)'),
    Poly(t**2 + 2*x*t + 1 + x**2, t, domain='ZZ[x]'), D, x, t, z) == \
        ([(Poly(z**2 + S(1)/4, z, domain='QQ'), Poly(t + x + 2*z, t,
        domain='ZZ[x,z]'))], True)


def test_integrate_hypertangent_polynomial():
    D = Poly(t**2 + 1, t)
    assert integrate_hypertangent_polynomial(Poly(t**2 + x*t + 1, t), D, x, t) == \
        (Poly(t, t, domain='ZZ'), Poly(x/2, t, domain='QQ[x]'))
    assert integrate_hypertangent_polynomial(Poly(t**5, t), Poly(a*(t**2 + 1), t), x, t) == \
        (Poly(1/(4*a)*t**4 - 1/(2*a)*t**2, t, domain='ZZ(a)'),
        Poly(1/(2*a), t, domain='ZZ(a)'))

def test_integrate_nonlinear_no_specials():
    a, d, = Poly(x**2*t**5 + x*t**4 - nu**2*t**3 - x*(x**2 + 1)*t**2 -(x**2 -
    nu**2)*t - x**5/4, t), Poly(x**2*t**4 + x**2*(x**2 + 2)*t**2 + x**2 +x**4 + x**6/4, t)
    D = Poly(-t**2 - t/x - (1 - nu**2/x**2), t)
    # f(x) == phi_nu(x), the logarithmic derivative of J_v, the Bessel function,
    # which has no specials (see Chapter 5, note 4 of Bronstein's book).
    f = Function('phi_nu')
    assert integrate_nonlinear_no_specials(a, d, D, x, t, f) == \
        (-log(1 + f(x)**2 + x**2/2)/2 - (4 + x**2)/(4 + 2*x**2 + 4*f(x)**2), True)
    assert integrate_nonlinear_no_specials(Poly(t, t), Poly(1, t), D, x, t, f) == \
        (0, False)
