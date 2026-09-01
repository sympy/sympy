from __future__ import annotations

from sympy.polys.domains.finitefield import GF
from sympy.polys.rings import ring
from sympy.polys.domains.integerring import ZZ
from sympy.polys.sparsegcd import smp_cofactors_FF, smp_cofactors_ZZ


def test_smp_cofactors_ZZ():
    # This gcd computation gets dispatched to zippel's algorithm.
    R, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9 = ring("x0, x1, x2, x3, x4, x5, x6, x7, x8, x9", ZZ)
    h = (x0**6 + x1**6 + x2**6 + x3**6 + x4**6 + x5**6 + x6**6 + x7**6 + x8**6 + x9**6)
    a = 1 + x0**6 + 2*x1**6 + 3*x2**6 + 4*x3**6 + 5*x4**6
    b = 1 + 7*x5**6 + 8*x6**6 + 9*x7**6 + 10*x8**6 + 11*x9**6

    f = h*a
    g = h*b

    gcd, cff, cfg = smp_cofactors_ZZ(f, g, 10)

    assert (R(gcd), R(cff), R(cfg)) == (h, a, b)


def test_smp_cofactors_FF():
    # This gcd computation gets dispatched to the modular zippel gcd algorithm.
    K = GF(100019)
    R, x, y, z = ring("x, y, z", K)
    h = x**3*y*z + y**2 + x**2 +5
    a = (x - 1)
    b = (x + 1)
    f = h*a
    g = h*b
    gcd, cff, cfg = smp_cofactors_FF(f, g, 3, K)

    assert (R(gcd), R(cff), R(cfg)) == (h, a, b)
