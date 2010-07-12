"""Most of these tests come from the examples in Bronstein's book."""
from sympy import Poly, S, Function, log, symbols, exp, tan, Integral
from sympy.integrals.risch import (gcdex_diophantine, derivation, splitfactor,
    splitfactor_sqf, canonical_representation, hermite_reduce,
    polynomial_reduce, residue_reduce, integrate_hyperexponential,
    integrate_hypertangent_polynomial, integrate_nonlinear_no_specials)
from sympy.utilities.pytest import XFAIL, skip

from sympy.abc import x, t, nu, z, a
t1, t2 = symbols('t1 t2')

def test_gcdex_diophantine():
    assert gcdex_diophantine(Poly(x**4 - 2*x**3 - 6*x**2 + 12*x + 15),
    Poly(x**3 + x**2 - 4*x - 4), Poly(x**2 - 1)) == \
        (Poly((-x**2 + 4*x - 3)/5), Poly((x**3 - 7*x**2 + 16*x - 10)/5))

def test_derivation():
    p = Poly(4*x**4*t**5 + (-4*x**3 - 4*x**4)*t**4 + (-3*x**2 + 2*x**3)*t**3 +
        (2*x + 7*x**2 + 2*x**3)*t**2 + (1 - 4*x - 4*x**2)*t - 1 + 2*x, t)
    D = [Poly(1, x), Poly(-t**2 - 3/(2*x)*t + 1/(2*x), t)]
    assert derivation(p, D, [x, t]) == Poly(-20*x**4*t**6 + (2*x**3 + 16*x**4)*t**5 +
        (21*x**2 + 12*x**3)*t**4 + (7*x/2 - 25*x**2 - 12*x**3)*t**3 +
        (-5 - 15*x/2 + 7*x**2)*t**2 - (3 - 8*x - 10*x**2 - 4*x**3)/(2*x)*t +
        (1 - 4*x**2)/(2*x), t)
    assert derivation(Poly(1, t), D, [x, t]) == Poly(0, t)
    assert derivation(Poly(t, t), D, [x, t]) == D[1]
    assert derivation(Poly(t**2 + 1/x*t + (1 - 2*x)/(4*x**2), t), D, [x, t]) == \
        Poly(-2*t**3 - 4/x*t**2 - (5 - 2*x)/(2*x**2)*t - (1 - 2*x)/(2*x**3), t, domain='ZZ(x)')
    # TODO: Add tests for multiple extensions and coefficientD
    assert derivation(Poly(x*t*t1, t), [Poly(1, x), Poly(1/x, t1), Poly(t, t)], [x, t1, t]) == \
        Poly(t*t1 + x*t*t1 + t, t)
    assert derivation(Poly(x*t*t1, t), [Poly(1, x), Poly(1/x, t1), Poly(t, t)], [x, t1, t],
    coefficientD=True) ==  Poly((1 + t1)*t, t)
    assert derivation(Poly(x, x), [Poly(1, x)], [x]) == Poly(1, x)

def test_splitfactor():
    p = Poly(4*x**4*t**5 + (-4*x**3 - 4*x**4)*t**4 + (-3*x**2 + 2*x**3)*t**3 +
        (2*x + 7*x**2 + 2*x**3)*t**2 + (1 - 4*x - 4*x**2)*t - 1 + 2*x, t, field=True)
    D = [Poly(1, x), Poly(-t**2 - 3/(2*x)*t + 1/(2*x), t)]
    assert splitfactor(p, D, [x, t]) == (Poly(4*x**4*t**3 + (-8*x**3 - 4*x**4)*t**2 +
        (4*x**2 + 8*x**3)*t - 4*x**2, t), Poly(t**2 + 1/x*t + (1 - 2*x)/(4*x**2), t, domain='ZZ(x)'))
    assert splitfactor(Poly(x, t), D, [x, t]) == (Poly(x, t), Poly(1, t))
    r = Poly(-4*x**4*z**2 + 4*x**6*z**2 - z*x**3 - 4*x**5*z**3 + 4*x**3*z**3 + x**4 + z*x**5 - x**6, t)
    D = [Poly(1, x), Poly(1/x, t)]
    assert splitfactor(r, D, [x, t], coefficientD=True) == \
        (Poly(x*z - x**2 - z*x**3 + x**4, t), Poly(-x**2 + 4*x**2*z**2, t))
    assert splitfactor_sqf(r, D, [x, t], coefficientD=True) == \
        (((Poly(x*z - x**2 - z*x**3 + x**4, t), 1),), ((Poly(-x**2 + 4*x**2*z**2, t), 1),))
    assert splitfactor(Poly(0, t), D, [x, t]) == (Poly(0, t), Poly(1, t))
    assert splitfactor_sqf(Poly(0, t), D, [x, t]) == (((Poly(0, t), 1),), ())

def test_canonical_representation():
    D = [Poly(1, x), Poly(1 + t**2, t)]
    assert canonical_representation(Poly(x - t, t), Poly(t**2, t), D, [x, t]) == \
        (Poly(0, t), (Poly(0, t),
        Poly(1, t)), (Poly(-t + x, t),
        Poly(t**2, t)))
    D = [Poly(1, x), Poly(t**2 + 1, t)]
    assert canonical_representation(Poly(t**5 + t**3 + x**2*t + 1, t),
    Poly((t**2 + 1)**3, t), D, [x, t]) == \
        (Poly(0, t), (Poly(t**5 + t**3 + x**2*t + 1, t),
        Poly(t**6 + 3*t**4 + 3*t**2 + 1, t)), (Poly(0, t), Poly(1, t)))

def test_hermite_reduce():
    D = [Poly(1, x), Poly(t**2 + 1, t)]
    assert hermite_reduce(Poly(x - t, t), Poly(t**2, t), D, [x, t]) == \
        ((Poly(-x, t), Poly(t, t)), (Poly(0, t), Poly(1, t)), (Poly(-x, t), Poly(1, t)))
    D = [Poly(1, x), Poly(-t**2 - t/x - (1 - nu**2/x**2), t)]
    assert hermite_reduce(Poly(x**2*t**5 + x*t**4 - nu**2*t**3 - x*(x**2 + 1)*t**2 -
    (x**2 - nu**2)*t - x**5/4, t), Poly(x**2*t**4 + x**2*(x**2 + 2)*t**2 + x**2 +
    x**4 + x**6/4, t), D, [x, t]) == \
        ((Poly(-1 - x**2/4, t), Poly(t**2 + 1 + x**2/2, t)),
        (Poly((-2*nu**2 - x**4)/(2*x**2)*t - (1 + x**2)/x, t),
        Poly(t**2 + 1 + x**2/2, t)), (Poly(t + 1/x, t), Poly(1, t)))
    D = [Poly(1, x), Poly(1/x, t)]
    assert hermite_reduce(Poly(-t**2 + 2*t + 2, t),
    Poly(-x*t**2 + 2*x*t - x, t), D, [x, t]) == \
        ((Poly(3, t), Poly(t - 1, t)), (Poly(0, t), Poly(1, t)), (Poly(1, t), Poly(x, t)))
    assert hermite_reduce(Poly(-x**2*t**6 + (-1 - 2*x**3 + x**4)*t**3 +
    (-3 - 3*x**4)*t**2 - 2*x*t - x - 3*x**2, t),
    Poly(x**4*t**6 - 2*x**2*t**3 + 1, t), D, [x, t]) == \
        ((Poly(x**5*t + x**2 + x**6, t), Poly(x**5*t**3 - x**3, t)), (Poly(0, t),
        Poly(1, t)), (Poly(-1, t), Poly(x**2, t)))

def test_polynomial_reduce():
    D = [Poly(1, x), Poly(1 + t**2, t)]
    assert polynomial_reduce(Poly(1 + x*t + t**2, t), D, [x, t]) == \
        (Poly(t, t), Poly(x*t, t))
    assert polynomial_reduce(Poly(0, t), D, [x, t]) == \
        (Poly(0, t), Poly(0, t))

def test_residue_reduce():
    a = Poly(2*t**2 - t - x**2, t)
    d = Poly(t**3 - x**2*t, t)
    D = [Poly(1, x), Poly(1/x, t)]
    assert residue_reduce(a, d, D, [x, t], z, invert=False) == \
        ([(Poly(z**2 - S(1)/4, z), Poly((1 + 3*x*z - 6*z**2 -
        2*x**2 + 4*x**2*z**2)*t - x*z + x**2 + 2*x**2*z**2 - 2*z*x**3, t))], False)
    assert residue_reduce(a, d, D, [x, t], z, invert=True) == \
        ([(Poly(z**2 - S(1)/4, z), Poly(t + 2*x*z, t))], False)
    assert residue_reduce(Poly(-2/x, t), Poly(t**2 - 1, t,), D, [x, t], z, invert=False) == \
        ([(Poly(z**2 - 1, z), Poly(-z*t - 1, t))], True)
    assert residue_reduce(Poly(-2/x, t), Poly(t**2 - 1, t), D, [x, t], z, invert=True) == \
        ([(Poly(z**2 - 1, z), Poly(t + z, t))], True)
    D = [Poly(1, x), Poly(-t**2 - t/x - (1 - nu**2/x**2), t)]
    # TODO: Skip or make faster
    assert residue_reduce(Poly((-2*nu**2 - x**4)/(2*x**2)*t - (1 + x**2)/x, t),
    Poly(t**2 + 1 + x**2/2, t), D, [x, t], z) == \
        ([(Poly(1, z), Poly(t, t)), (Poly(z + S(1)/2, z), Poly(t**2 + 1 + x**2/2, t))], True)
    D = [Poly(1, x), Poly(1 + t**2, t)]
    assert residue_reduce(Poly(-2*x*t + 1 - x**2, t),
    Poly(t**2 + 2*x*t + 1 + x**2, t), D, [x, t], z) == \
        ([(Poly(z**2 + S(1)/4, z), Poly(t + x + 2*z, t))], True)

def test_integrate_hyperexponential():
    # TODO: Add tests for integrate_hyperexponential() from the book
    a = Poly((1 + 2*t1 + t1**2 + 2*t1**3)*t**2 + (1 + t1**2)*t + 1 + t1**2, t)
    d = Poly(1, t)
    D = [Poly(1, x), Poly(1 + t1**2, t1), Poly(t*(1 + t1**2), t)]
    assert integrate_hyperexponential(a, d, D, [x, t1, t], [lambda x: exp(tan(x)), tan]) == \
        (exp(2*tan(x))*tan(x) + Integral(1 + tan(x)**2, x) + exp(tan(x)), True)
        # (exp(2*tan(x))*tan(x) + tan(x) + exp(tan(x)), True)

    a = Poly(t, t)
    d = Poly(1, t)
    D = [Poly(1, x), Poly(2*x*t, t)]

    assert integrate_hyperexponential(a, d, D, [x, t], [lambda x: exp(x**2)]) == (0, False)

    D = [Poly(1, x), Poly(t, t)]
    assert integrate_hyperexponential(a, d, D, [x, t], [exp]) == (exp(x), True)

    a = Poly(25*t**6 - 10*t**5 + 7*t**4 - 8*t**3 + 13*t**2 + 2*t - 1, t)
    d = Poly(25*t**6 + 35*t**4 + 11*t**2 + 1, t)
    assert integrate_hyperexponential(a, d, D, [x, t], [exp]) == \
        (-(55 - 50*exp(x))/(25 + 125*exp(2*x)) + Integral(-1, x) + log(1 + exp(2*x)), True)
        # (-(55 - 50*exp(x))/(25 + 125*exp(2*x)) - x + log(1 + exp(2*x)), True)

def test_integrate_hypertangent_polynomial():
    D = [Poly(1, x), Poly(t**2 + 1, t)]
    assert integrate_hypertangent_polynomial(Poly(t**2 + x*t + 1, t), D, [x, t]) == \
        (Poly(t, t), Poly(x/2, t))
    D = [Poly(1, x), Poly(a*(t**2 + 1), t)]
    assert integrate_hypertangent_polynomial(Poly(t**5, t), D, [x, t]) == \
        (Poly(1/(4*a)*t**4 - 1/(2*a)*t**2, t), Poly(1/(2*a), t))

def test_integrate_nonlinear_no_specials():
    a, d, = Poly(x**2*t**5 + x*t**4 - nu**2*t**3 - x*(x**2 + 1)*t**2 -(x**2 -
    nu**2)*t - x**5/4, t), Poly(x**2*t**4 + x**2*(x**2 + 2)*t**2 + x**2 +x**4 + x**6/4, t)
    D = [Poly(1, x), Poly(-t**2 - t/x - (1 - nu**2/x**2), t)]
    # f(x) == phi_nu(x), the logarithmic derivative of J_v, the Bessel function,
    # which has no specials (see Chapter 5, note 4 of Bronstein's book).
    f = Function('phi_nu')
    assert integrate_nonlinear_no_specials(a, d, D, [x, t], [f]) == \
        (-log(1 + f(x)**2 + x**2/2)/2 - (4 + x**2)/(4 + 2*x**2 + 4*f(x)**2), True)
    assert integrate_nonlinear_no_specials(Poly(t, t), Poly(1, t), D, [x, t], [x, f]) == \
        (0, False)
