from __future__ import annotations
from sympy.external.gmpy import MPZ
from sympy.matrices.dense import Matrix
from sympy.polys.densebasic import dup_to_dict
from sympy.polys.rings import ring
from sympy.polys.domains import ZZ
from sympy.polys.zippel import (
    from_newt_to_poly, incremental_newton_interp, lag_basis, skeleton_sorter, smp_LC_wrt_last, smp_chinese_remainder_reconstruction_multivariate,
    smp_deg_wrt_last, smp_gf_gcd, smp_primitive_wrt_last,
    smp_trivial_gcd, vandermonde_interp)


def test_smp_gf_gcd():
    R, x = ring("x", ZZ)
    p = MPZ(7)
    f = (x**2 + 2).trunc_ground(p)
    g = (x**2 + 3).trunc_ground(p)
    assert smp_gf_gcd(dict(f), dict(g), p, R.domain) == dict(R.one)

    p = MPZ(11)
    f = (x**2 + 5*x + 4).trunc_ground(p)
    g = (x**2 + 6*x + 8).trunc_ground(p)
    expected = (x + 4).trunc_ground(p)
    assert smp_gf_gcd(dict(f), dict(g), p, R.domain) == dict(expected)

    # gcd in Z[x] would be not monic
    p = MPZ(13)
    f = (2*x**2 + 5*x + 3).trunc_ground(p)
    g = (2*x**2 - x - 6).trunc_ground(p)
    expected = (x - 5).trunc_ground(p)
    assert smp_gf_gcd(dict(f), dict(g), p, R.domain) == dict(expected)


def test_smp_trivial_gcd():
    R, x, y, z = ring("x, y, z", ZZ)
    n = R.ngens
    dom = R.domain

    assert smp_trivial_gcd(dict(R.zero), dict(R.zero), n, dom) == ({}, {}, {})

    g = 2*x*y*z + x**2
    assert smp_trivial_gcd(dict(R.zero), dict(g), n, dom) == (
        dict(g), dict(R.zero), dict(R.one))

    g_neg = -3*x**3*y - z
    assert smp_trivial_gcd(dict(R.zero), dict(g_neg), n, dom) == (
        dict(-g_neg), dict(R.zero), dict(-R.one))

    f = 5*y**2 + x
    assert smp_trivial_gcd(dict(f), dict(R.zero), n, dom) == (
        dict(f), dict(R.one), dict(R.zero))

    f_const = R.ground_new(24)
    g_const = R.ground_new(36)
    result = smp_trivial_gcd(dict(f_const), dict(g_const), n, dom)
    assert result == (
        dict(R.ground_new(12)),
        dict(R.ground_new(2)),
        dict(R.ground_new(3)))

    assert smp_trivial_gcd(dict(x**2 + y), dict(x + z), n, dom) is None


def test_smp_primitive_wrt_last():
    R, x, y, z = ring("x, y, z", ZZ)
    p = MPZ(7)
    C = (3*z**2 + 2*z + 1).trunc_ground(p)
    P = (x**2*y - x*z + 2).trunc_ground(p)
    f = (C * P).trunc_ground(p)
    C_expected = (5 * C).trunc_ground(p)
    P_expected = (3 * P).trunc_ground(p)
    contf, ppf = smp_primitive_wrt_last(dict(f), R.ngens, R.domain, p)

    expected_contf = {
        (monom[-1],): R.domain.convert(coeff % p)
        for monom, coeff in dict(C_expected).items()
    }
    assert dup_to_dict(contf, R.domain) == expected_contf
    assert ppf == dict(P_expected)
    assert (R.from_dict(ppf) * C_expected).trunc_ground(p) == f


def test_smp_LC_wrt_last():
    R, x, y = ring("x, y", ZZ)
    f = x**2 * y**2 + 3 * x**2 + 5 * x * y + 7
    R_, y_ = ring("y", ZZ)
    expected_lc = y_**2 + 3

    assert smp_LC_wrt_last(dict(f), R.ngens, R.domain) == dict(expected_lc)


def test_smp_deg_wrt_last():
    R, x, y, z = ring("x, y, z", ZZ)
    f1 = x**3 * y**2 * z**5 + x**4 * y * z + x**2 * y**3 * z**8
    assert smp_deg_wrt_last(dict(f1), R.ngens) == (4, 1)

    f2 = x**2 * y * z**3 + x * y**100 * z**5
    assert smp_deg_wrt_last(dict(f2), R.ngens) == (2, 1)


def test_smp_chinese_remainder_reconstruction_multivariate():
    R, x, y = ring("x, y", ZZ)
    p, q = MPZ(3), MPZ(5)

    hp = x**3*y - x**2 - 1
    hq = -x**3*y - 2*x*y**2 + 2

    hpq_dict = smp_chinese_remainder_reconstruction_multivariate(
        dict(hp), dict(hq), p, q, R.domain, R.ngens)
    hpq = R.from_dict(hpq_dict)

    assert hpq.trunc_ground(p) == hp
    assert hpq.trunc_ground(q) == hq

    R, x, y, z = ring("x, y, z", ZZ)
    p, q = MPZ(6), MPZ(5)

    hp = 3*x**4 - y**3*z + z
    hq = -2*x**4 + z

    hpq_dict = smp_chinese_remainder_reconstruction_multivariate(
        dict(hp), dict(hq), p, q, R.domain, R.ngens)
    hpq = R.from_dict(hpq_dict)

    assert hpq.trunc_ground(p) == hp
    assert hpq.trunc_ground(q) == hq


def test_vandermonde_interpolate():
    p = ZZ(100003)
    k = ZZ.map([3, 6, 12, 33])
    A = Matrix([[el**j for j in range(len(k))] for el in k])
    x_list = [12, 2, 1, 27]
    x = Matrix(x_list)
    v = [ZZ(int(el)) for el in A * x]
    v_t = [ZZ(int(el)) for el in A.T * x]

    bas = lag_basis(k, p)
    sol_t = vandermonde_interp(bas, v_t, p)
    sol = vandermonde_interp(bas, v, p, trans=False)

    assert sol == sol_t == x_list

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

    assert from_newt_to_poly(x, v, p) == [3, 99, 17, 96]


def test_skeleton_sorter():
    R, x, y, z, w = ring("x, y, z, w", ZZ)
    G = 4 * x**3 * y**2 * w - 2 * x**3 * z**5 * w**2 + 7 * x * y * z * w
    S, h, monic, pseudomonic = skeleton_sorter(dict(G))

    assert S == {
        1: [((1, 1, 1, 1), [(0, 1), (1, 1), (2, 1)])],
        3: [((3, 2, 0, 1), [(0, 2), (2, 1)]), ((3, 0, 5, 2), [(1, 5), (2, 2)])],
    }
    assert h == [[7], [4], [-2]]
    assert monic == pseudomonic == False

    R, x, y, z = ring("x, y, z", ZZ)
    G = 5 * x**2 * y**2 + 7 * x * y**5 * z**3 + 8 * x * z**4
    S, h, monic, pseudomonic = skeleton_sorter(G)

    assert S == {
        2: [((2, 2, 0), [(0, 2)])],
        1: [((1, 5, 3), [(0, 5), (1, 3)]), ((1, 0, 4), [(1, 4)])],
    }
    assert list(S.keys()) == [2, 1]
    assert h == [[5], [7], [8]]
    assert monic == pseudomonic == True
