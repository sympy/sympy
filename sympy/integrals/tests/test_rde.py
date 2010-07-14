"""Most of these tests come from the examples in Bronstein's book."""
from sympy import Poly, S, symbols, oo
from sympy.integrals.risch import NonElementaryIntegral
from sympy.integrals.rde import (order_at, weak_normalizer, normal_denom,
    bound_degree, spde, solve_poly_rde, no_cancel_equal, rischDE)
from sympy.utilities.pytest import raises
from sympy.abc import x, t, z, n

t0, t1, t2 = symbols('t0, t1, t2')

def test_order_at():
    a = Poly(t**4, t)
    b = Poly((t**2 + 1)**3*t, t)
    p1 = Poly(t, t)
    p2 = Poly(1 + t**2, t)
    assert order_at(a, p1, t) == 4
    assert order_at(b, p1, t) == 1
    assert order_at(a, p2, t) == 0
    assert order_at(b, p2, t) == 3
    raises(ValueError, "order_at(Poly(0, t), Poly(t, t), t)")

def test_weak_normalizer():
    a = Poly((1 + x)*t**5 + 4*t**4 + (-1 - 3*x)*t**3 - 4*t**2 + (-2 + 2*x)*t, t)
    d = Poly(t**4 - 3*t**2 + 2, t)
    D = [Poly(1, x), Poly(t, t)]
    r = weak_normalizer(a, d, D, [x, t])
    assert r == (Poly(t**5 - t**4 - 4*t**3 + 4*t**2 + 4*t - 4, t),
        (Poly((1 + x)*t**2 + x*t, t), Poly(t + 1, t)))
    assert weak_normalizer(r[1][0], r[1][1], D, [x, t]) == (Poly(1, t), r[1])
    r = weak_normalizer(Poly(1 + t**2), Poly(t**2 - 1, t), D, [x, t], z)
    assert r == (Poly(t**4 - 2*t**2 + 1, t), (Poly(-3*t**2 + 1, t), Poly(t**2 - 1, t)))
    assert weak_normalizer(r[1][0], r[1][1], D, [x, t]) == (Poly(1, t), r[1])
    D = [Poly(1, x), Poly(1 + t**2)]
    r = weak_normalizer(Poly(1 + t**2), Poly(t, t), D, [x, t], z)
    assert r == (Poly(t, t), (Poly(0, t), Poly(1, t)))
    assert weak_normalizer(r[1][0], r[1][1], D, [x, t]) == (Poly(1, t), r[1])

def test_normal_denom():
    raises(NonElementaryIntegral, """normal_denom(Poly(1, x), Poly(1, x),
    Poly(1, x), Poly(x, x), [Poly(1, x)], [x])""")
    fa, fd = Poly(t**2 + 1, t), Poly(1, t)
    ga, gd = Poly(1, t), Poly(t**2, t)
    D = [Poly(1, x), Poly(t**2 + 1, t)]
    assert normal_denom(fa, fd, ga, gd, D, [x, t]) == \
        (Poly(t, t), (Poly(t**3 - t**2 + t - 1, t), Poly(1, t)), (Poly(1, t),
        Poly(1, t)), Poly(t, t))

def test_special_denom():
    # TODO: add more tests here
    D = [Poly(t, t)]
    assert special_denom(Poly(1, t), Poly(t**2, t), Poly(1, t), Poly(t**2 - 1, t),
    Poly(t, t), D, [t]) == \
        (Poly(1, t), Poly(t**2 - 1, t), Poly(t**2 - 1, t), Poly(t, t))

def test_bound_degree():
    # Primitive (TODO)

    # Base (TODO)
    D = [Poly(1, x)]
    assert bound_degree(Poly(1, x), Poly(-2*x, x), Poly(1, x), D, [x]) == 0

    # Exp (TODO)

    # Nonlinear
    D = [Poly(1, x), Poly(t**2 + 1, t)]
    assert bound_degree(Poly(t, t), Poly((t - 1)*(t**2 + 1), t), Poly(1, t),
        D, [x, t]) == 0

def test_spde():
    raises(NonElementaryIntegral, "spde(Poly(t, t), Poly((t - 1)*(t**2 + 1), " +
        "t), Poly(1, t), [Poly(1, x), Poly(1 + t**2, t)], 0, [x, t])")
    D = [Poly(1, x), Poly(t, t)]
    assert spde(Poly(t**2 + x*t*2 + x**2, t), Poly(t**2/x**2 + (2/x - 1)*t, t),
    Poly(t**2/x**2 + (2/x - 1)*t, t), D, 0, [x, t]) == \
        (Poly(0, t), Poly(0, t), 0, Poly(0, t), Poly(1, t))
    D = [Poly(1, x), Poly(t0/x**2, t0), Poly(1/x, t)]
    assert spde(Poly(t**2, t), Poly(-t**2/x**2 - 1/x, t),
    Poly((2*x - 1)*t**4 + (t0 + x)/x*t**3 - (t0 + 4*x**2)/(2*x)*t**2 + x*t, t),
    D, 3, [x, t0, t]) == \
        (Poly(0, t), Poly(0, t), 0, Poly(0, t),
        Poly(t0*t**2/2 + x**2*t**2 - x**2*t, t))
    D = [Poly(1, x)]
    assert spde(Poly(x**2 + x + 1, x), Poly(-2*x - 1, x), Poly(x**5/2 +
    3*x**4/4 + x**3 - x**2 + 1, x), D, 4, [x]) == \
        (Poly(0, x), Poly(x/2 - S(1)/4, x), 2, Poly(x**2 + x + 1, x), Poly(5*x/4, x))
    assert spde(Poly(x**2 + x + 1, x), Poly(-2*x - 1, x), Poly(x**5/2 +
    3*x**4/4 + x**3 - x**2 + 1, x), D, n, [x]) == \
        (Poly(0, x), Poly(x/2 - S(1)/4, x), -2 + n, Poly(x**2 + x + 1, x), Poly(5*x/4, x))

def test_solve_poly_rde_no_cancel():
    # deg(b) large
    assert solve_poly_rde(Poly(t**2 + 1, t), Poly(t**3 + (x + 1)*t**2 + t + x + 2, t),
    [Poly(1, x), Poly(1 + t**2, t)], oo, [x, t]) == Poly(t + x, t)
    # deg(b) small
    assert solve_poly_rde(Poly(0, x), Poly(x/2 - S(1)/4, x), [Poly(1, x)], oo, [x]) == \
        Poly(x**2/4 - x/4, x)
    D = [Poly(1, x), Poly(t**2 + 1, t)]
    assert solve_poly_rde(Poly(2, t), Poly(t**2 + 2*t + 3, t), D, 1, [x, t]) == \
        Poly(t + 1, t, x)
    # deg(b) == deg(D) - 1
    D = [Poly(1, x), Poly(t**2 + 1, t)]
    assert no_cancel_equal(Poly(1 - t, t),
    Poly(t**3 + t**2 - 2*x*t - 2*x, t), D, oo, [x, t]) == \
        (Poly(t**2, t), 1, Poly((-2 - 2*x)*t - 2*x, t))

def test_rischDE():
    # TODO: Add more tests for rischDE, including ones from the text
    assert rischDE(Poly(-2*x, x), Poly(1, x), Poly(1 - 2*x - 2*x**2, x),
    Poly(1, x), [Poly(1, x)], [x]) == \
        (Poly(x + 1, x), Poly(1, x))
