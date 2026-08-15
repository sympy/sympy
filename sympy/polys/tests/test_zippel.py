from __future__ import annotations
from unittest.mock import patch
from sympy.external.gmpy import MPZ
from sympy.matrices.dense import Matrix
from sympy.polys.densebasic import dup_to_dict
from sympy.polys.rings import ring
from sympy.polys.domains import ZZ
from sympy.polys.zippel import (
    from_newt_to_poly, incremental_newton_interp, lag_basis, skeleton_sorter, smp_LC_wrt_last, smp_chinese_remainder_reconstruction_multivariate,
    smp_deg_wrt_last, smp_gf_gcd, smp_primitive_wrt_last,
    smp_trivial_gcd, smp_zippel_gcd, smp_zippel_gcd_mod, smp_zippel_interp_mod, vandermonde_interp)


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
    f = 5*y**2 - x
    assert smp_trivial_gcd(dict(f), {}, n, dom)

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


def test_smp_zippel_gcd():
    # trivial tests
    R, x, y = ring("x, y", ZZ)
    f = (x + 1)*(x + y)
    g = {}
    gcd, cff, cfg = smp_zippel_gcd(f, g, 2)
    assert (R(gcd), R(cff), R(cfg)) == (f, R.one, R.zero)

    R, x, y = ring("x, y", ZZ)
    # monic tests
    f = (x + 1)*(x + y)
    g = (x - 1)*(x + y)
    gcd, cff, cfg = smp_zippel_gcd(f, g, 2)
    assert (R(gcd), R(cff), R(cfg)) == (x + y, x+1,x-1)

    R, x, y, z = ring("x, y, z", ZZ)
    f = (x + 1)*(x + 3*y**2*z + 2*z)
    g = (x - 1)*(x + 3*y**2*z + 2*z)
    gcd, cff, cfg = smp_zippel_gcd(f, g, 3)
    assert (R(gcd), R(cff), R(cfg)) == (x + 3*y**2*z + 2*z, x+1,x-1)

    # pseudomonic tests
    R, x, y = ring("x, y", ZZ)
    f = (x + 1)*(x*y + y)
    g = (x - 1)*(x*y + y)
    gcd, cff, cfg = smp_zippel_gcd(f, g, 2)
    assert (R(gcd), R(cff), R(cfg)) == (x*y + y, x+1,x-1)

    R, x, y, z = ring("x, y, z", ZZ)
    f = (x + 1)*(x*y**2*z + 3*y**2*z + 2*z)
    g = (x - 1)*(x*y**2*z + 3*y**2*z + 2*z)
    gcd, cff, cfg = smp_zippel_gcd(f, g, 3)
    assert (R(gcd), R(cff), R(cfg)) == (x*y**2*z + 3*y**2*z + 2*z, x+1,x-1)

    # not monic tests
    R, x, y = ring("x, y", ZZ)
    f = (x**2+ x + 1)*(x*y + y)
    g = (x**2+ 2*x + 1)*(x*y + y)
    gcd, cff, cfg = smp_zippel_gcd(f, g, 2)
    assert (R(gcd), R(cff), R(cfg)) == (x*y + y, x**2+ x + 1, x**2+ 2*x + 1)

    R, x, y, z = ring("x, y, z", ZZ)
    f = (x + 1)*(x + 2*x*y + 3*y**2*z + 2*z)
    g = (x - 1)*(x + 2*x*y + 3*y**2*z + 2*z)
    gcd, cff, cfg = smp_zippel_gcd(f, g, 3)
    assert (R(gcd), R(cff), R(cfg)) == (x + 2*x*y + 3*y**2*z + 2*z, x+1,x-1)

    f = (x**2+1)*(2*x**2*y + 3*x**2*z + 13*x + 7*y*z)
    g = (x**2 -1)*(2*x**2*y + 3*x**2*z + 13*x + 7*y*z)
    gcd, cff, cfg = smp_zippel_gcd(f, g, 3)
    assert (R(gcd), R(cff), R(cfg)) == (2*x**2*y + 3*x**2*z + 13*x + 7*y*z, x**2+1,x**2-1)

    R, x, y, z, w = ring("x, y, z, w", ZZ)
    h = x**2*(z**2 + z + y) + x*(y**2 + y + 2*z) + (y*z +2*y +3*z)
    f = (x + w)*h
    g = (x - w)*h
    gcd, cff, cfg = smp_zippel_gcd(f, g, 4)
    assert (R(gcd), R(cff), R(cfg)) == (h, x+w,x-w)

    # test with first prime canceling out leading coeff
    R, x, y = ring("x, y", ZZ)
    f = (1000000007*x + 1)*(x +y)
    g = (1000000007*x - 1)*(x+y)
    gcd, cff, cfg = smp_zippel_gcd(f, g, 2)
    assert (R(gcd), R(cff), R(cfg)) == (x +y, 1000000007*x + 1,1000000007*x - 1)

    # test unlucky first prime
    R, x, y = ring("x, y", ZZ)
    f = (x + y)*(x + 1)
    g = (x + y)*(x + 1000000007 + 1)
    gcd, cff, cfg = smp_zippel_gcd(f, g, 2)
    assert (R(gcd), R(cff), R(cfg)) == (x + y, x + 1, x + 1000000007 + 1)

    # test unlucky second prime
    R, x, y = ring("x, y", ZZ)
    f = (x + y)*(x + 1)
    g = (x + y)*(x + 1000000009 + 1)
    gcd, cff, cfg = smp_zippel_gcd(f, g, 2)
    assert (R(gcd), R(cff), R(cfg)) == (x + y, x + 1, x + 1000000009 + 1)

    # test unlucky sequence of evaluations in smp_zippel_gcd_mod
    R, x, y = ring("x, y", ZZ)
    q = (y - 1)*(y - 2)*(y - 3)
    f = x*q + 1
    g = x*q + y
    with patch("sympy.polys.zippel.randint", side_effect=[1, 2, 3, 4, 4],):
        gcd, cff, cfg = smp_zippel_gcd(
            dict(f), dict(g), 2)
    assert (R(gcd), R(cff), R(cfg)) == (R.one, f, g)

    # unlucky first evaluation leads to wrong skeleton, that leads to inconsistent system
    R, x, y, z = ring("x, y, z", ZZ)
    H2 = 1 + y + y**2 + y**4
    H1 = 1 + 2*y + 3*y**2 + 4*y**4
    H0 = 1 + 3*y + 5*y**2 + 7*y**4
    h = (
        x**2*(H2 + (z - 1)*y**3)
        + x*(H1 + (z - 1)*y**5)
        + H0 + (z - 1)*y**6
    )
    f = (x + 1)*h
    g = (x + 2)*h
    with patch(
        "sympy.polys.zippel._random_eval_point",
        side_effect=[1, 17, 23, 2, 3, 4, 5, 6, 7, 8]):
        gcd, cff, cfg = smp_zippel_gcd(f, g, 3)
    assert (R(gcd), R(cff), R(cfg)) == (h, x + 1, x + 2)

    # returned gcd's LC with positive sign
    R, x, y = ring("x, y", ZZ)
    f = (x + 1)*(-x**2 + y)
    g = (x + y)*(-x**2 + y)
    gcd, cff, cfg = smp_zippel_gcd(f, g, 2)
    assert (R(gcd), R(cff), R(cfg)) == (x**2 - y, - x - 1, -x - y)

    R, x, y = ring("x, y", ZZ)
    p1 = ZZ(1000000007)
    p2 = ZZ(1000000009)
    h = x + y
    f = ((p1*p2 - 1)*x + 1)*h
    g = ((p1*p2 - 1)*x - 1)*h
    gcd, cff, cfg = smp_zippel_gcd(f, g, 2)
    assert (R(gcd), R(cff), R(cfg)) == (h , (p1*p2 - 1)*x + 1, (p1*p2 - 1)*x - 1)

    # the first 2 primes are unlucky and lead to failing trial division
    R, x, y, z = ring("x, y, z", ZZ)
    p1 = ZZ(1000000007)
    p2 = ZZ(1000000009)
    f = (x + 1)*(x**2 + x*(y*p1*p2 + z) + y**2 + z**2)
    g = (x + y)*(x**2 + x*(y*p1*p2 + z) + y**2 + z**2)
    gcd, cff, cfg = smp_zippel_gcd(f, g, 3)
    assert (R(gcd), R(cff), R(cfg)) == (x**2 + x*(y*p1*p2 + z) + y**2 + z**2, x + 1, x + y)

    # the gcd has a power of x with zero coefficient
    R, x, y, z = ring("x, y, z", ZZ)
    f = (x + 1)*(x**2 + x*z)
    g = (x + y)*(x**2 + x*z)
    gcd, cff, cfg = smp_zippel_gcd(f, g, 3)
    assert (R(gcd), R(cff), R(cfg)) == (x**2 + x*z, x + 1, x + y)

def test_smp_zippel_gcd_mod():
    # prime choice always leads to LC canceling out after evaluation of y
    R, x, y = ring("x, y", ZZ)
    f = x*(y + 1) + 1
    g = x*(y + 1) + y
    assert(smp_zippel_gcd_mod(f, g, ZZ(2), 2) is None)

    # coprime inputs lead to early exit
    gcd = smp_zippel_gcd_mod(x+y, x-y, ZZ(43), 2)
    assert R(gcd) == 1

    # wrong gcd is determined by recursive call due to unlucky evaluation of z
    R, x, y, z = ring("x, y, z", ZZ)
    f = (x + y - z)*(x**2*y + x*(z-3)+ y*(z-2)+(z-5))
    g = (x + y + z)*(x**2*y + x*(z-3)+ y*(z-2)+(z-5))
    with patch("sympy.polys.zippel._random_eval_point", side_effect=[2, 3, 5],):
        gcd = smp_zippel_gcd_mod(f, g, ZZ(10009), 3)
    assert gcd is None

    # smp_zippel_gcd_mod returns None in recursive call
    R, x, y, z, w = ring("x, y, z, w", ZZ)
    f = (x + y - z*w)*(x**2*y + x*(z-3)+ y*(z-2)+(z-5) + w)
    g = (x + y + z*w)*(x**2*y + x*(z-3)+ y*(z-2)+(z-5) + w)
    with patch("sympy.polys.zippel._random_eval_point", side_effect=[2, 3, 5, 2, 3, 5, 2, 3, 5],):
        gcd = smp_zippel_gcd_mod(f, g, ZZ(10009), 4)
    assert (gcd is None)

    # y*(z - 2) incorrectly interpolated as 2*(z - 2) by repeatedly
    # evaluating y in 2, leads to fail final division test
    R, x, y, z = ring("x, y, z", ZZ)
    h = x**2*y + y*(z - 2) + 1
    f = (x + y)*h
    g = (x - y)*h
    with (
        patch(
            "sympy.polys.zippel._random_eval_point",
            side_effect=[2, 1, 2, 3, 3, 4, 5, 6, 7, 8],
        ),
        patch(
            "sympy.polys.zippel._random_eval_tuple",
            return_value=(2,),
        ),
    ):
        gcd = smp_zippel_gcd_mod(f, g, ZZ(10009), 3)
    assert gcd is None

    # test unlucky first evaluation
    R, x, y = ring("x, y", ZZ)
    f = x + 17
    g = x + y
    with (
        patch(
            "sympy.polys.zippel._random_eval_point",
            side_effect=[17, 3, 4, 5, 6, 7, 8],
        ),
    ):
        gcd = smp_zippel_gcd_mod(f, g, ZZ(10009), 2)
    assert gcd is None

    f = (x + 1)*(x**2 + (y - 17)*x)
    g = (x - 1)*(x**2 + (y - 17)*x)
    with (
        patch(
            "sympy.polys.zippel._random_eval_point",
            side_effect=[17, 3, 4, 5, 6, 7, 8],
        ),
    ):
        gcd = smp_zippel_gcd_mod(f, g, ZZ(10009), 2)
    assert gcd is None

    # test unlucky first evaluation
    R, x, y, z = ring("x, y, z", ZZ)
    f = x + y + z
    g = x + y + 17
    with (
        patch(
            "sympy.polys.zippel._random_eval_point",
            side_effect=[17, 3, 4, 5, 6, 7, 8],
        ),
    ):
        gcd = smp_zippel_gcd_mod(f, g, ZZ(10009), 3)
    assert gcd is None

    # test repeated choice of same evaluation point
    R, x, y = ring("x, y", ZZ)
    f = (x + 1)*(x + y)
    g = (x - 1)*(x + y)
    with (
        patch(
            "sympy.polys.zippel.randint",
            side_effect=[3, 3, 3, 5, 6, 7, 8, 9],
        ),
    ):
        gcd = smp_zippel_gcd_mod(f, g, ZZ(10009), 2)
    assert R(gcd) == x + y

    # not-monic case with content fails, preprocessing needs to extract it
    R, x, y, z = ring("x, y, z", ZZ)
    h = x**2*(z**2 + y) + x*(z**2 + y) + (z**2 + y)
    f = (z - 1)*h
    g = (z + 1)*h
    gcd = smp_zippel_gcd_mod(f, g, ZZ(10009), 3)
    assert gcd is None

    # unlucky first evaluation leads to wrong skeleton, that leads to inconsistent system
    R, x, y, z = ring("x, y, z", ZZ)
    H2 = 1 + y + y**2 + y**4
    H1 = 1 + 2*y + 3*y**2 + 4*y**4
    H0 = 1 + 3*y + 5*y**2 + 7*y**4
    h = (
        x**2*(H2 + (z - 1)*y**3)
        + x*(H1 + (z - 1)*y**5)
        + H0 + (z - 1)*y**6
    )
    f = (x + 1)*h
    g = (x + 2)*h
    with patch(
        "sympy.polys.zippel._random_eval_point",
        side_effect=[1, 17, 23]):
        gcd = smp_zippel_gcd_mod(f, g, ZZ(10009), 3)
    assert gcd is None

    # unlucky evaluation leads to content wrt x in non monic case
    R, x, y, z = ring("x, y, z", ZZ)
    f = (y + 1)*(x**2*(y + z) + x*(z -1))
    g = (y - 1)*(x**2*(y + z) + x*(z -1))
    with patch(
        "sympy.polys.zippel._random_eval_point",
        side_effect=[1, 17, 23]):
        gcd = smp_zippel_gcd_mod(f, g, ZZ(10009), 3)
    assert gcd is None

    # unlucky evaluation leads to content wrt x in monic case
    R, x, y, z = ring("x, y, z", ZZ)
    f = (y + 1)*(x**2*(y) + x*(z -1))
    g = (y - 1)*(x**2*(y) + x*(z -1))
    with patch(
        "sympy.polys.zippel._random_eval_point",
        side_effect=[1, 17, 23]):
        gcd = smp_zippel_gcd_mod(f, g, ZZ(10009), 3)
    assert gcd is None

def test_smp_zippel_interp_mod():
    # testing bad evaluation tuples
    R, x, y, z = ring("x, y, z", ZZ)
    h = x*y + x*z + 1
    f = (x + 1)*h
    g = (x + 2)*h
    G, _, monic, pseudomonic = skeleton_sorter(h)
    with patch(
        "sympy.polys.zippel._random_eval_tuple",
        side_effect=[(1, 100), (1, 1), (1, 2)]):
        result = smp_zippel_interp_mod(
            f,
            g,
            G,
            ZZ(101),
            monic,
            pseudomonic,
            3,
        )
    assert result == {(0, 0, 0): 34, (1, 1, 0): 34, (1, 0, 1): 34}

    # evaluation tuple becomes bad when an additional image is added
    R, x, y, z = ring("x, y, z", ZZ)
    P = y + z
    Q = y**2 + y*z + z**2
    h = (x + 1)*P + x**2*Q
    f = (x + 1)*h
    g = (x - 1)*h
    G, _, monic, pseudomonic = skeleton_sorter(h)
    with patch("sympy.polys.zippel._random_eval_tuple",
        side_effect=[(ZZ(1), ZZ(2)), (ZZ(1), ZZ(5))]):
        gcd = smp_zippel_interp_mod(f, g, G, ZZ(13), monic, pseudomonic, 3)
    assert R(gcd) == R(gcd).LC*h

    # evaluation tuple cancels the leading coefficient of B but not A
    R, x, y = ring("x, y", ZZ)
    h = x + 1
    f = (x + 1)*(x + 2)
    g = (x + 1)*((y - 2)*x + 1)
    G, _, monic, pseudomonic = skeleton_sorter(h)
    with patch(
        "sympy.polys.zippel._random_eval_tuple",
        side_effect=[(ZZ(2),), (ZZ(3),)]):
        gcd = smp_zippel_interp_mod(f, g, G, ZZ(101), monic, pseudomonic, 2)
    assert R(gcd) == h

    # here the first number of equations (z) chosen "a priori" won't be enough, and will be updated
    R, x, y, z = ring("x, y, z", ZZ)
    p = ZZ(101)
    P = y + z
    Q = y**2 + y*z + z**2
    h = (x + 1)*P + x**2*Q
    f = (x + 1)*h
    g = (x - 1)*h
    G, _, monic, pseudomonic = skeleton_sorter(h)
    gcd = smp_zippel_interp_mod(f, g, G, p, monic, pseudomonic, 3)
    assert R(gcd)== R(gcd).LC*h
