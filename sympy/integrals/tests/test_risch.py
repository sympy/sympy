"""Most of these tests come from the examples in Bronstein's book."""
from sympy import (Poly, S, Function, log, symbols, exp, tan, Integral, sqrt,
    Symbol, Lambda, sin)
from sympy.integrals.risch import (gcdex_diophantine, frac_in, as_poly_1t,
    derivation, splitfactor, splitfactor_sqf, canonical_representation,
    hermite_reduce, polynomial_reduce, residue_reduce, residue_reduce_to_basic,
    integrate_primitive, integrate_hyperexponential_polynomial,
    integrate_hyperexponential, integrate_hypertangent_polynomial,
    integrate_nonlinear_no_specials, integer_powers, DifferentialExtension,
    risch_integrate)
from sympy.utilities.pytest import XFAIL, skip, raises

from sympy.abc import x, t, nu, z, a, y
t0, t1, t2 = symbols('t0, t1, t2')

def test_gcdex_diophantine():
    assert gcdex_diophantine(Poly(x**4 - 2*x**3 - 6*x**2 + 12*x + 15),
    Poly(x**3 + x**2 - 4*x - 4), Poly(x**2 - 1)) == \
        (Poly((-x**2 + 4*x - 3)/5), Poly((x**3 - 7*x**2 + 16*x - 10)/5))

def test_frac_in():
    assert frac_in(Poly((x + 1)/x*t, t), x) == \
        (Poly(t*x + t, x), Poly(x, x))
    assert frac_in((x + 1)/x*t, x) == \
        (Poly(t*x + t, x), Poly(x, x))
    assert frac_in((Poly((x + 1)/x*t, t), Poly(t + 1, t)), x) == \
        (Poly(t*x + t, x), Poly((1 + t)*x, x))
    raises(ValueError, "frac_in((x + 1)/log(x)*t, x)")
    assert frac_in(Poly((2 + 2*x + x*(1 + x))/(1 + x)**2, t), x, cancel=True) == \
        (Poly(x + 2, x), Poly(x + 1, x))

def test_as_poly_1t():
    assert as_poly_1t(2/t + t, t, z) in [
        Poly(t + 2*z, t, z), Poly(t + 2*z, z, t)]
    assert as_poly_1t(2/t + 3/t**2, t, z) in [
        Poly(2*z + 3*z**2, t, z), Poly(2*z + 3*z**2, z, t)]
    assert as_poly_1t(2/((exp(2) + 1)*t), t, z) in [
        Poly(2/(exp(2) + 1)*z, t, z), Poly(2/(exp(2) + 1)*z, z, t)]
    assert as_poly_1t(2/((exp(2) + 1)*t) + t, t, z) in [
        Poly(t + 2/(exp(2) + 1)*z, t, z), Poly(t + 2/(exp(2) + 1)*z, z, t)]

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
    # Test basic option
    assert derivation((x + 1)/(x - 1), [Poly(1, x)], [x], basic=True) == -2/(1 - 2*x + x**2)
    assert derivation((t + 1)/(t - 1), [Poly(t, t)], [t], basic=True) == -2*t/(1 - 2*t + t**2)
    assert derivation(t + 1, [Poly(t, t)], [t], basic=True) == t

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
    # TODO: Fix this domain='EX' bug
    assert hermite_reduce(Poly(x**2*t**5 + x*t**4 - nu**2*t**3 - x*(x**2 + 1)*t**2 -
    (x**2 - nu**2)*t - x**5/4, t), Poly(x**2*t**4 + x**2*(x**2 + 2)*t**2 + x**2 +
    x**4 + x**6/4, t), D, [x, t]) == \
        ((Poly(-1 - x**2/4, t, domain='EX'), Poly(t**2 + 1 + x**2/2, t, domain='EX')),
        (Poly((2*nu**2 + x**4)/-(2*x**2)*t - (1 + x**2)/x, t, domain='EX', expand=False),
        Poly(t**2 + 1 + x**2/2, t, domain='EX')), (Poly(t + 1/x, t, domain='EX'),
        Poly(1, t, domain='EX')))
    D = [Poly(1, x), Poly(1/x, t)]
    assert hermite_reduce(Poly(-t**2 + 2*t + 2, t),
    Poly(-x*t**2 + 2*x*t - x, t), D, [x, t]) == \
        ((Poly(3, t), Poly(t - 1, t)), (Poly(0, t), Poly(1, t)), (Poly(1, t), Poly(x, t)))
    assert hermite_reduce(Poly(-x**2*t**6 + (-1 - 2*x**3 + x**4)*t**3 +
    (-3 - 3*x**4)*t**2 - 2*x*t - x - 3*x**2, t),
    Poly(x**4*t**6 - 2*x**2*t**3 + 1, t), D, [x, t]) == \
        ((Poly(x**5*t + x**2 + x**6, t), Poly(x**5*t**3 - x**3, t)), (Poly(0, t),
        Poly(1, t)), (Poly(-1, t), Poly(x**2, t)))
    assert hermite_reduce(Poly((-2 + 3*x)*t0**3 + (-1 + x)*t0**2 +
    (-4*x + 2*x**2)*t0 + x**2, t0), Poly(x*t0**6 - 4*x**2*t0**5 +
    6*x**3*t0**4 - 4*x**4*t0**3 + x**5*t0**2, t0), D, [x, t0]) == \
        ((Poly(t0**2 + t0/3 + x, t0), Poly(t0**4 - 3*x*t0**3 + 3*x**2*t0**2 -
        x**3*t0, t0)), (Poly(0, t0), Poly(1, t0)), (Poly(0, t0), Poly(1, t0)))

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
    ans = residue_reduce(Poly(-2/x, t), Poly(t**2 - 1, t), D, [x, t], z, invert=True)
    assert ans == ([(Poly(z**2 - 1, z), Poly(t + z, t))], True)
    assert residue_reduce_to_basic(ans[0], [x, t], z, [log]) == -log(-1 + log(x)) + log(1 + log(x))

    D = [Poly(1, x), Poly(-t**2 - t/x - (1 - nu**2/x**2), t)]
    # TODO: Skip or make faster
    assert residue_reduce(Poly((-2*nu**2 - x**4)/(2*x**2)*t - (1 + x**2)/x, t),
    Poly(t**2 + 1 + x**2/2, t), D, [x, t], z) == \
        ([(Poly(z + S(1)/2, z, domain='QQ'), Poly(t**2 + 1 + x**2/2, t, domain='EX'))], True)
    D = [Poly(1, x), Poly(1 + t**2, t)]
    assert residue_reduce(Poly(-2*x*t + 1 - x**2, t),
    Poly(t**2 + 2*x*t + 1 + x**2, t), D, [x, t], z) == \
        ([(Poly(z**2 + S(1)/4, z), Poly(t + x + 2*z, t))], True)
    D = [Poly(1, x), Poly(t, t)]
    assert residue_reduce(Poly(t, t), Poly(t + sqrt(2), t), D, [x, t], z) == \
        ([(Poly(z - 1, z), Poly(t + sqrt(2), t))], True)

def test_integrate_hyperexponential():
    # TODO: Add tests for integrate_hyperexponential() from the book
    a = Poly((1 + 2*t1 + t1**2 + 2*t1**3)*t**2 + (1 + t1**2)*t + 1 + t1**2, t)
    d = Poly(1, t)
    D = [Poly(1, x), Poly(1 + t1**2, t1), Poly(t*(1 + t1**2), t)]
    assert integrate_hyperexponential(a, d, D, [x, t1, t], [lambda x: exp(tan(x)), tan]) == \
        (exp(2*tan(x))*tan(x) + exp(tan(x)), 1 + t1**2, True)
        # exp(2*tan(x))*tan(x) + tan(x) + exp(tan(x))
    a = Poly((t1**3 + (x + 1)*t1**2 + t1 + x + 2)*t, t)
    assert integrate_hyperexponential(a, d, D, [x, t1, t], [lambda x: exp(tan(x)), tan]) == \
        ((x + tan(x))*exp(tan(x)), 0, True)

    a = Poly(t, t)
    d = Poly(1, t)
    D = [Poly(1, x), Poly(2*x*t, t)]

    assert integrate_hyperexponential(a, d, D, [x, t], [lambda x: exp(x**2)]) == \
        (0, Integral(exp(x**2), x), False)

    D = [Poly(1, x), Poly(t, t)]
    assert integrate_hyperexponential(a, d, D, [x, t], [exp]) == \
        (exp(x), 0, True)

    a = Poly(25*t**6 - 10*t**5 + 7*t**4 - 8*t**3 + 13*t**2 + 2*t - 1, t)
    d = Poly(25*t**6 + 35*t**4 + 11*t**2 + 1, t)
    assert integrate_hyperexponential(a, d, D, [x, t], [exp]) == \
        (-(55 - 50*exp(x))/(25 + 125*exp(2*x)) + log(1 + exp(2*x)), -1, True)
        # -(55 - 50*exp(x))/(25 + 125*exp(2*x)) - x + log(1 + exp(2*x))
    D = [Poly(1, x), Poly(t0, t0), Poly(t0*t, t)]
    assert integrate_hyperexponential(Poly(2*t0*t**2, t), Poly(1, t), D, [x, t0, t],
    [lambda x: exp(exp(x)), exp]) == \
        (exp(2*exp(x)), 0, True)

    D = [Poly(1, x), Poly(t0, t0), Poly(-t0*t, t)]
    assert integrate_hyperexponential(Poly(-27*exp(9) - 162*t0*exp(9) +
    27*x*t0*exp(9), t), Poly((36*exp(18) + x**2*exp(18) - 12*x*exp(18))*t, t),
    D, [x, t0, t], [lambda x: exp(exp(x)), exp]) == \
        (27*exp(-exp(x))/(-6*exp(9) + x*exp(9)), 0, True)

    D = [Poly(1, x), Poly(t, t)]
    assert integrate_hyperexponential(Poly(x**2/2*t, t), Poly(1, t), D, [x, t], [exp]) == \
        ((2 - 2*x + x**2)*exp(x)/2, 0, True)
    assert integrate_hyperexponential(Poly(1 + t, t), Poly(t, t), D, [x, t], [exp]) == \
        (- exp(-x), 1, True) # x - exp(-x)
    assert integrate_hyperexponential(Poly(x, t), Poly(t + 1, t), D, [x, t], [exp]) == \
        (0, Integral(x/(1 + exp(x)), x), False)

def test_integrate_hyperexponential_polynomial():
    # Without proper cancelation within integrate_hyperexponential_polynomial(),
    # this will take a long time to complete, and will return a complicated
    # expression
    p = Poly((-28*x**11*t0 - 6*x**8*t0 + 6*x**9*t0 - 15*x**8*t0**2 +
        15*x**7*t0**2 + 84*x**10*t0**2 - 140*x**9*t0**3 - 20*x**6*t0**3 +
        20*x**7*t0**3 - 15*x**6*t0**4 + 15*x**5*t0**4 + 140*x**8*t0**4 -
        84*x**7*t0**5 - 6*x**4*t0**5 + 6*x**5*t0**5 + x**3*t0**6 - x**4*t0**6 +
        28*x**6*t0**6 - 4*x**5*t0**7 + x**9 - x**10 + 4*x**12)/(-8*x**11*t0 +
        28*x**10*t0**2 - 56*x**9*t0**3 + 70*x**8*t0**4 - 56*x**7*t0**5 +
        28*x**6*t0**6 - 8*x**5*t0**7 + x**4*t0**8 + x**12)*t1**2 +
        (-28*x**11*t0 - 12*x**8*t0 + 12*x**9*t0 - 30*x**8*t0**2 +
        30*x**7*t0**2 + 84*x**10*t0**2 - 140*x**9*t0**3 - 40*x**6*t0**3 +
        40*x**7*t0**3 - 30*x**6*t0**4 + 30*x**5*t0**4 + 140*x**8*t0**4 -
        84*x**7*t0**5 - 12*x**4*t0**5 + 12*x**5*t0**5 - 2*x**4*t0**6 +
        2*x**3*t0**6 + 28*x**6*t0**6 - 4*x**5*t0**7 + 2*x**9 - 2*x**10 +
        4*x**12)/(-8*x**11*t0 + 28*x**10*t0**2 - 56*x**9*t0**3 +
        70*x**8*t0**4 - 56*x**7*t0**5 + 28*x**6*t0**6 - 8*x**5*t0**7 +
        x**4*t0**8 + x**12)*t1 + (-2*x**2*t0 + 2*x**3*t0 + x*t0**2 -
        x**2*t0**2 + x**3 - x**4)/(-4*x**5*t0 + 6*x**4*t0**2 - 4*x**3*t0**3 +
        x**2*t0**4 + x**6), t1, z, expand=False)
    D = [Poly(1, x), Poly(1/x, t0), Poly(2*x*t1, t1)]
    T = [x, t0, t1]
    assert integrate_hyperexponential_polynomial(p, D, T, z) == (
        Poly((x - t0)*t1**2 + (-2*t0 + 2*x)*t1, t1), Poly(-2*x*t0 + x**2 +
        t0**2, t1), True)

def test_integrate_primitive():
    D = [Poly(1, x), Poly(1/x, t)]
    assert integrate_primitive(Poly(t, t), Poly(1, t), D, [x, t], [log]) == \
        (x*log(x), -1, True) # (x*log(x) - x, True)
    assert integrate_primitive(Poly(x, t), Poly(t, t), D, [x, t], [log]) == \
        (0, Integral(x/log(x), x), False)

    D = [Poly(1, x), Poly(1/x, t1), Poly(1/(x + 1), t2)]
    assert integrate_primitive(Poly(t1, t2), Poly(t2, t2), D, [x, t1, t2],
    [lambda x: log(x + 1), log]) == \
        (0, Integral(log(x)/log(1 + x), x), False)

    D = [Poly(1, x), Poly(1/x, t1), Poly(1/(x*t1), t2)]
    assert integrate_primitive(Poly(t2, t2), Poly(t1, t2), D, [x, t1, t2],
    [lambda x: log(log(x)), log]) == \
        (0, Integral(log(log(x))/log(x), x), False)

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

def test_integer_powers():
    assert integer_powers([x, x/2, x**2 + 1, 2*x/3], index=True) == \
        [(x/6, [(0, 6), (1, 3), (3, 4)]), (1 + x**2, [(2, 1)])]
    assert integer_powers([x, x/2, x**2 + 1, 2*x/3], index=False) == \
        [(x/6, [(x, 6), (x/2, 3), (2*x/3, 4)]), (1 + x**2, [(1 + x**2, 1)])]

def test_DifferentialExtension():
    # XXX: These might be different with different .arg orderings
    i = Symbol('i')

    # Exponentials
    assert DifferentialExtension(exp(x) + exp(x**2), x,
    dummy=False)._important_attrs == \
        (Poly(t1 + t0, t1), Poly(1, t1), [Poly(1, x,), Poly(t0, t0),
        Poly(2*x*t1, t1)], [x, t0, t1], [Lambda(i, exp(i**2)),
        Lambda(i, exp(i))], [], [1, 2], [x, x**2], [], [])
    assert DifferentialExtension(exp(x) + exp(2*x), x,
    dummy=False)._important_attrs == \
        (Poly(t0**2 + t0, t0), Poly(1, t0), [Poly(1, x), Poly(t0, t0)], [x, t0],
        [Lambda(i, exp(i))], [], [1], [x], [], [])
    assert DifferentialExtension(exp(x) + exp(x/2), x,
    dummy=False)._important_attrs == \
        (Poly(t0**2 + t0, t0), Poly(1, t0), [Poly(1, x), Poly(t0/2, t0)],
        [x, t0], [Lambda(i, exp(i/2))], [], [1], [x/2], [], [])
    assert DifferentialExtension(exp(x) + exp(x**2) + exp(x + x**2), x,
    dummy=False)._important_attrs == \
        (Poly((1 + t0)*t1 + t0, t1), Poly(1, t1), [Poly(1, x), Poly(t0, t0),
        Poly(2*x*t1, t1)], [x, t0, t1], [Lambda(i, exp(i**2)),
        Lambda(i, exp(i))], [], [1, 2], [x, x**2], [], [])
    assert DifferentialExtension(exp(x) + exp(x**2) + exp(x + x**2 + 1), x,
    dummy=False)._important_attrs == \
        (Poly((1 + S.Exp1*t0)*t1 + t0, t1), Poly(1, t1), [Poly(1, x),
        Poly(t0, t0), Poly(2*x*t1, t1)], [x, t0, t1], [Lambda(i, exp(i**2)),
        Lambda(i, exp(i))], [], [1, 2], [x, x**2], [], [])
    assert DifferentialExtension(exp(x) + exp(x**2) + exp(x/2 + x**2), x,
    dummy=False)._important_attrs == \
        (Poly(t1**2 + t0*t1 + t0, t1), Poly(1, t1), [Poly(1, x),
        Poly(2*x*t0, t0), Poly(t1/2, t1)], [x, t0, t1], [Lambda(i, exp(i/2)),
        Lambda(i, exp(i**2))], [], [1, 2], [x**2, x/2], [], [])
    assert DifferentialExtension(exp(x) + exp(x**2) + exp(x/2 + x**2 + 3), x,
    dummy=False)._important_attrs == \
        (Poly(t1**2 + t0*exp(3)*t1 + t0, t1), Poly(1, t1), [Poly(1, x),
        Poly(2*x*t0, t0), Poly(t1/2, t1)], [x, t0, t1], [Lambda(i, exp(i/2)),
        Lambda(i, exp(i**2))], [], [1, 2], [x**2, x/2], [], [])

    # Logarithms
    assert DifferentialExtension(log(x)*log(x + 1)*log(2*x**2 + 2*x), x,
    dummy=False)._important_attrs == \
        (Poly(t0*t1**2 + (-t0*log(2) - t0**2)*t1, t1), Poly(1, t1),
        [Poly(1, x), Poly(1/(1 + x), t0),
        Poly((1 + 2*x)/(x + x**2), t1, expand=False)], [x, t0, t1],
        [Lambda(i, log(2*i + 2*i**2)), Lambda(i, log(1 + i))], [], [], [],
        [1, 2], [x + 1, 2*x**2 + 2*x])
    assert DifferentialExtension(x**x*log(x), x,
    dummy=False)._important_attrs == \
        (Poly(t0*t1, t1), Poly(1, t1), [Poly(1, x), Poly(1/x, t0),
        Poly((1 + t0)*t1, t1)], [x, t0, t1], [Lambda(i, exp(t0*i)),
        Lambda(i, log(i))], [(exp(x*log(x)), x**x)], [2], [t0*x], [1], [x])

    # Test handle_first
    assert DifferentialExtension(exp(x)*log(x), x,
    handle_first='log', dummy=False)._important_attrs == \
        (Poly(t0*t1, t1), Poly(1, t1), [Poly(1, x), Poly(1/x, t0),
        Poly(t1, t1)], [x, t0, t1], [Lambda(i, exp(i)), Lambda(i, log(i))],
        [], [2], [x], [1], [x])
    assert DifferentialExtension(exp(x)*log(x), x,
    handle_first='exp', dummy=False)._important_attrs == \
        (Poly(t0*t1, t1), Poly(1, t1), [Poly(1, x), Poly(t0, t0),
        Poly(1/x, t1)], [x, t0, t1], [Lambda(i, log(i)), Lambda(i, exp(i))],
        [], [1], [x], [2], [x])

    # This one must have the log first, regardless of what we set it to
    # (because the log is inside of the exponential: x**x == exp(x*log(x)))
    assert DifferentialExtension(-x**x*log(x)**2 + x**x - x**x/x, x,
    handle_first='exp', dummy=False)._important_attrs == \
    DifferentialExtension(-x**x*log(x)**2 + x**x - x**x/x, x,
    handle_first='log', dummy=False)._important_attrs == \
        (Poly((-1 + x - x*t0**2)*t1, t1), Poly(x, t1),
        [Poly(1, x), Poly(1/x, t0), Poly((1 + t0)*t1, t1)], [x, t0, t1],
        [Lambda(i, exp(t0*i)), Lambda(i, log(i))], [(exp(x*log(x)), x**x)],
        [2], [t0*x], [1], [x])

    # Test 'unimportant' attributes
    DE = DifferentialExtension(exp(x)*log(x), x, dummy=False, handle_first='exp')
    assert DE.f == exp(x)*log(x)
    assert DE.newf == t0*t1
    assert DE.x == x
    assert DE.cases == ['primitive', 'exp', 'base']

    assert DE.level == -1
    assert DE.t == t1 == DE.T[DE.level]
    assert DE.d == Poly(1/x, t1) == DE.D[DE.level]
    raises(ValueError, 'DE.increment_level()')
    DE.decrement_level()
    assert DE.level == -2
    assert DE.t == t0 == DE.T[DE.level]
    assert DE.d == Poly(t0, t0) == DE.D[DE.level]
    DE.decrement_level()
    assert DE.level == -3
    assert DE.t == x == DE.T[DE.level] == DE.x
    assert DE.d == Poly(1, x) == DE.D[DE.level]
    raises(ValueError, 'DE.decrement_level()')
    DE.increment_level()
    DE.increment_level()
    assert DE.level == -1
    assert DE.t == t1 == DE.T[DE.level]
    assert DE.d == Poly(1/x, t1) == DE.D[DE.level]

    # Test the extension flag
    raises(ValueError, "DifferentialExtension(None, None, extension={'T':[x, t]})")
    DE = DifferentialExtension(None, None, extension={'D':[Poly(1, x), Poly(t, t)]})
    assert DE._important_attrs == (None, None, [Poly(1, x), Poly(t, t)], [x, t],
        None, None, None, None, None, None)
    assert DE.d == Poly(t, t)
    assert DE.t == t
    assert DE.level == -1
    assert DE.cases == ['exp', 'base']

    DE = DifferentialExtension(None, None, extension={'D':[Poly(1, x),
        Poly(t, t)], 'E_K':[1], 'E_args':[x], 'L_K':[], 'L_args':[]})
    assert DE._important_attrs == (None, None, [Poly(1, x), Poly(t, t)], [x, t],
        None, None, [1], [x], [], [])

    # Odd ends
    assert DifferentialExtension(sin(y)*exp(x), x,
    dummy=False)._important_attrs == \
        (Poly(sin(y)*t0, t0, domain='ZZ[sin(y)]'), Poly(1, t0, domain='ZZ'),
        [Poly(1, x, domain='ZZ'), Poly(t0, t0, domain='ZZ')], [x, t0],
        [Lambda(i, exp(i))], [], [1], [x], [], [])
    raises(NotImplementedError, "DifferentialExtension(sin(x), x)")

    # Rothstein's integral
    f = (2581284541*exp(x) + 1757211400)/(39916800*exp(3*x) +
    119750400*exp(x)**2 + 119750400*exp(x) + 39916800)*exp(1/(exp(x) + 1) - 10*x)
    assert DifferentialExtension(f, x, dummy=False)._important_attrs == \
        (Poly((1757211400 + 2581284541*t0)*t1, t1), Poly(39916800 +
        119750400*t0 + 119750400*t0**2 + 39916800*t0**3, t1),
        [Poly(1, x), Poly(t0, t0), Poly(-(10 + 21*t0 + 10*t0**2)/(1 + 2*t0 +
        t0**2)*t1, t1, domain='ZZ(t0)')], [x, t0, t1],
        [Lambda(i, exp(-10*i + 1/(1 + t0))), Lambda(i, exp(i))], [], [1, 2],
        [x, -10*x + 1/(1 + t0)], [], [])

def test_risch_integrate():
    assert risch_integrate(t0*exp(x), x) == t0*exp(x)

    # From my GSoC writeup
    assert risch_integrate((1 + 2*x**2 + x**4 + 2*x**3*exp(2*x**2))/
    (x**4*exp(x**2) + 2*x**2*exp(x**2) + exp(x**2)), x) == \
        Integral(exp(-x**2), x) + exp(x**2)/(1 + x**2)
