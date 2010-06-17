"""Most of these tests come from the examples in Bronstein's book."""
from sympy import Poly
from sympy.integrals.rde import (weak_normalizer, )
from sympy.abc import x, t, z

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
