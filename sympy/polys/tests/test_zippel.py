from __future__ import annotations
from sympy.polys.rings import ring
from sympy.polys.domains import ZZ
from sympy.polys.zippel import (_LC, _chinese_remainder_reconstruction_multivariate, _deg,
    _gf_gcd, _trivial_gcd, _primitive, skeleton_sorter, from_newt_to_poly,
    incremental_newton_interp)


def test_gf_gcd():
    R, x = ring("x", ZZ)
    p = 7
    f = (x**2 + 2).trunc_ground(p)
    g = (x**2 + 3).trunc_ground(p)
    assert _gf_gcd(f, g, p) == R.one

    p = 11
    f = (x**2 + 5*x + 4).trunc_ground(p)
    g = (x**2 + 6*x + 8).trunc_ground(p)
    assert _gf_gcd(f, g, p) == (x + 4).trunc_ground(p)

    # gcd in Z[x] would be not monic
    p = 13
    f = (2*x**2 + 5*x + 3).trunc_ground(p)
    g = (2*x**2 - x - 6).trunc_ground(p)
    assert _gf_gcd(f, g, p) == (x - 5).trunc_ground(p)


def test_trivial_gcd():
    R, x, y, z = ring("x, y, z", ZZ)

    assert _trivial_gcd(R.zero, R.zero) == (R.zero, R.zero, R.zero)

    g = 2*x*y*z + x**2
    assert _trivial_gcd(R.zero, g) == (g, R.zero, R.one)

    g_neg = -3*x**3*y - z
    assert _trivial_gcd(R.zero, g_neg) == (-g_neg, R.zero, -R.one)

    f = 5*y**2 + x
    assert _trivial_gcd(f, R.zero) == (f, R.one, R.zero)

    f_const = R.ground_new(24)
    g_const = R.ground_new(36)
    h, cff, cfg = _trivial_gcd(f_const, g_const)
    assert h == R.ground_new(12)
    assert cff == R.ground_new(2)
    assert cfg == R.ground_new(3)

    assert _trivial_gcd(x**2 + y, x + z) is None


def test_primitive():
    R, x, y, z = ring("x, y, z", ZZ)
    p = 7
    C = (3*z**2 + 2*z + 1).trunc_ground(p)
    P = (x**2*y - x*z + 2).trunc_ground(p)
    f = (C * P).trunc_ground(p)
    C_expected = (5 * C).trunc_ground(p)
    P_expected = (3 * P).trunc_ground(p)
    contf, ppf = _primitive(f, p)

    assert contf.set_ring(R) == C_expected
    assert ppf == P_expected
    assert (contf.set_ring(R) * ppf).trunc_ground(p) == f


def test_LC():
    R, x, y = ring("x, y", ZZ)
    f = x**2 * y**2 + 3 * x**2 + 5 * x * y + 7
    R_, y_ = ring("y", ZZ)
    expected_lc = y_**2 + 3

    assert _LC(f) == expected_lc


def test_deg():
    R, x, y, z = ring("x, y, z", ZZ)
    f1 = x**3 * y**2 * z**5 + x**4 * y * z + x**2 * y**3 * z**8
    assert _deg(f1) == (4, 1)

    f2 = x**2 * y * z**3 + x * y**100 * z**5
    assert _deg(f2) == (2, 1)


def test_chinese_remainder_reconstruction_multivariate():
    R, x, y = ring("x, y", ZZ)
    p, q = 3, 5

    hp = x**3*y - x**2 - 1
    hq = -x**3*y - 2*x*y**2 + 2

    hpq = _chinese_remainder_reconstruction_multivariate(hp, hq, p, q)

    assert hpq.trunc_ground(p) == hp
    assert hpq.trunc_ground(q) == hq

    T, z = ring("z", R)
    p, q = 3, 7

    hp = (x*y + 1)*z**2 + x
    hq = (x**2 - 3*y)*z + 2

    hpq = _chinese_remainder_reconstruction_multivariate(hp, hq, p, q)

    assert hpq.trunc_ground(p) == hp
    assert hpq.trunc_ground(q) == hq


def test_incremental_newton_interp():
    # Polynomial to interpolate: P(x) = 3x^3 - 2x^2 + 17x - 5
    x = ZZ.map([0, 1, 2])
    v = ZZ.map([96, 18, 7])
    xk = ZZ(3)
    uk = ZZ(8)
    p = ZZ(101)

    assert incremental_newton_interp(x, v, xk, uk, p) == 3


def test_from_newt_to_poly():
    x = ZZ.map([0, 1, 2, 3])
    v = ZZ.map([96, 18, 7, 3])
    p = ZZ(101)

    assert from_newt_to_poly(x, v, p) == [96, 17, 99, 3]


def test_skeleton_sorter():
    R, x, y, z, w = ring("x, y, z, w", ZZ)
    G = 4*x**3*y**2*w - 2*x**3*z**5*w**2 + 7*x*y*z*w
    S, h, monic, pseudomonic = skeleton_sorter(dict(G))

    assert S == {
        1: [[(1, 1, 1, 1), (0, 1), (1, 1), (2, 1)]],
        3: [
            [(3, 2, 0, 1), (0, 2), (2, 1)],
            [(3, 0, 5, 2), (1, 5), (2, 2)]
        ]
    }
    assert h == [[7], [4], [-2]]
    assert monic == pseudomonic == False

    R, x, y, z = ring("x, y, z", ZZ)
    G = 5*x**2*y**2 + 7*x*y**5*z**3 + 8*x*z**4
    S, h, monic, pseudomonic = skeleton_sorter(G)

    assert S == {
        2: [[(2, 2, 0), (0, 2)]],
        1: [[(1, 5, 3), (0, 5), (1, 3)], [(1, 0, 4), (1, 4)]]
    }
    assert list(S.keys()) == [2, 1]
    assert h == [[5], [7], [8]]
    assert monic == pseudomonic == True
