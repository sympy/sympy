"""Most of these tests come from the examples in Bronstein's book."""
from sympy import Poly, S
from sympy.integrals.risch import NonElementaryIntegral
from sympy.integrals.rde import (order_at, weak_normalizer, normal_denominator,
    bound_degree, spde)
from sympy.utilities.pytest import raises
from sympy.abc import x, t, z, n

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
    a = Poly((1 + x)*t**5 + 4*t**4 + (-1 - 3*x)*t**3 - 4*t**2 + (-2 + 2*x)*t, t, domain='ZZ[x]')
    d = Poly(t**4 - 3*t**2 + 2, t, domain='ZZ')
    D = Poly(t, t)
    r = weak_normalizer(a, d, D, x, t)
    assert r == (Poly(t**5 - t**4 - 4*t**3 + 4*t**2 + 4*t - 4, t, domain='ZZ[x]'),
        (Poly((1 + x)*t**2 + x*t, t, domain='ZZ[x]'), Poly(t + 1, t, domain='ZZ[x]')))
    assert weak_normalizer(r[1][0], r[1][1], D, x, t) == (Poly(1, t), r[1])
    r = weak_normalizer(Poly(1 + t**2), Poly(t**2 - 1, t), D, x, t, z)
    assert r == (Poly(t**4 - 2*t**2 + 1, t, domain='ZZ'),
        (Poly(-3*t**2 + 1, t, domain='ZZ'), Poly(t**2 - 1, t, domain='ZZ')))
    assert weak_normalizer(r[1][0], r[1][1], D, x, t) == (Poly(1, t), r[1])
    D = Poly(1 + t**2)
    r = weak_normalizer(Poly(1 + t**2), Poly(t, t), D, x, t, z)
    assert r == (Poly(t, t, domain='ZZ'), (Poly(0, t, domain='ZZ'), Poly(1, t, domain='ZZ')))
    assert weak_normalizer(r[1][0], r[1][1], D, x, t) == (Poly(1, t), r[1])

def test_normal_denominator():
    raises(NonElementaryIntegral, """normal_denominator(Poly(1, t), Poly(1, t),
    Poly(1, t), Poly(t, t), Poly(1, t), x, t)""")
    fa, fd = Poly(t**2 + 1, t), Poly(1, t)
    ga, gd = Poly(1, t), Poly(t**2, t)
    D = Poly(t**2 + 1, t)
    assert normal_denominator(fa, fd, ga, gd, D, x, t) == \
        (Poly(t, t, domain='ZZ'), (Poly(t**3 - t**2 + t - 1, t, domain='ZZ'),
        Poly(1, t, domain='ZZ')), (Poly(1, t, domain='ZZ'), Poly(1, t, domain='ZZ')),
        Poly(t, t, domain='ZZ'))

def test_bound_degree():
    # Primitive (TODO)

    # Base (TODO)
    assert bound_degree(Poly(1, t), Poly(-2*t, t), Poly(1, t), Poly(1, t), x, t) == 0

    # Exp (TODO)

    # Nonlinear
    assert bound_degree(Poly(t, t), Poly((t - 1)*(t**2 + 1), t), Poly(1, t),
        Poly(t**2 + 1, t), x, t) == 0

def test_spde():
    raises(NonElementaryIntegral, "spde(Poly(t, t), Poly((t - 1)*(t**2 + 1), " +
        "t), Poly(1, t), Poly(1 + t**2, t), 0, x, t)")
    assert spde(Poly(t**2 + x*t*2 + x**2, t), Poly(t**2/x**2 + (2/x - 1)*t, t),
    Poly(t**2/x**2 + (2/x - 1)*t, t), Poly(t, t), 0, x, t) == \
        (Poly(0, t, domain='ZZ'), Poly(0, t, domain='ZZ'), 0,
        Poly(0, t, domain='ZZ(x)'), Poly(1, t, domain='ZZ(x)'))
    # TODO: add example 6.4.3, pg. 204-5 (requires support for multiple extensions)
    assert spde(Poly(t**2 + t + 1, t), Poly(-2*t - 1, t), Poly(t**5/2 +
    3*t**4/4 + t**3 - t**2 + 1, t), Poly(1, t), 4, x, t) == \
        (Poly(0, t, domain='ZZ'), Poly(t/2 - S(1)/4, t, domain='QQ'), 2,
        Poly(t**2 + t + 1, t, domain='ZZ'), Poly(5*t/4, t, domain='QQ'))
    assert spde(Poly(t**2 + t + 1, t), Poly(-2*t - 1, t), Poly(t**5/2 +
    3*t**4/4 + t**3 - t**2 + 1, t), Poly(1, t), n, x, t) == \
        (Poly(0, t, domain='ZZ'), Poly(t/2 - S(1)/4, t, domain='QQ'), -2 + n,
        Poly(t**2 + t + 1, t, domain='ZZ'), Poly(5*t/4, t, domain='QQ'))
