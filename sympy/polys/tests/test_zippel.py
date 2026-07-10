from __future__ import annotations
from sympy.external.gmpy import MPZ
from sympy.polys.densebasic import dup_to_dict
from sympy.polys.rings import ring
from sympy.polys.domains import ZZ
from sympy.polys.zippel import (
    dict_LC_wrt_last, dict_chinese_remainder_reconstruction_multivariate,
    dict_deg_wrt_last, dict_gf_gcd, dict_primitive_wrt_last,
    dict_trivial_gcd)


def test_dict_gf_gcd():
    R, x = ring("x", ZZ)
    p = MPZ(7)
    f = (x**2 + 2).trunc_ground(p)
    g = (x**2 + 3).trunc_ground(p)
    assert dict_gf_gcd(dict(f), dict(g), p, R.domain) == dict(R.one)

    p = MPZ(11)
    f = (x**2 + 5*x + 4).trunc_ground(p)
    g = (x**2 + 6*x + 8).trunc_ground(p)
    expected = (x + 4).trunc_ground(p)
    assert dict_gf_gcd(dict(f), dict(g), p, R.domain) == dict(expected)

    # gcd in Z[x] would be not monic
    p = MPZ(13)
    f = (2*x**2 + 5*x + 3).trunc_ground(p)
    g = (2*x**2 - x - 6).trunc_ground(p)
    expected = (x - 5).trunc_ground(p)
    assert dict_gf_gcd(dict(f), dict(g), p, R.domain) == dict(expected)


def test_dict_trivial_gcd():
    R, x, y, z = ring("x, y, z", ZZ)
    n = R.ngens
    dom = R.domain

    assert dict_trivial_gcd(dict(R.zero), dict(R.zero), n, dom) == ({}, {}, {})

    g = 2*x*y*z + x**2
    assert dict_trivial_gcd(dict(R.zero), dict(g), n, dom) == (
        dict(g), dict(R.zero), dict(R.one))

    g_neg = -3*x**3*y - z
    assert dict_trivial_gcd(dict(R.zero), dict(g_neg), n, dom) == (
        dict(-g_neg), dict(R.zero), dict(-R.one))

    f = 5*y**2 + x
    assert dict_trivial_gcd(dict(f), dict(R.zero), n, dom) == (
        dict(f), dict(R.one), dict(R.zero))

    f_const = R.ground_new(24)
    g_const = R.ground_new(36)
    result = dict_trivial_gcd(dict(f_const), dict(g_const), n, dom)
    assert result == (
        dict(R.ground_new(12)),
        dict(R.ground_new(2)),
        dict(R.ground_new(3)))

    assert dict_trivial_gcd(dict(x**2 + y), dict(x + z), n, dom) is None


def test_dict_primitive_wrt_last():
    R, x, y, z = ring("x, y, z", ZZ)
    p = MPZ(7)
    C = (3*z**2 + 2*z + 1).trunc_ground(p)
    P = (x**2*y - x*z + 2).trunc_ground(p)
    f = (C * P).trunc_ground(p)
    C_expected = (5 * C).trunc_ground(p)
    P_expected = (3 * P).trunc_ground(p)
    contf, ppf = dict_primitive_wrt_last(dict(f), R.ngens, R.domain, p)

    expected_contf = {
        (monom[-1],): R.domain.convert(coeff % p)
        for monom, coeff in dict(C_expected).items()
    }
    assert dup_to_dict(contf, R.domain) == expected_contf
    assert ppf == dict(P_expected)
    assert (R.from_dict(ppf) * C_expected).trunc_ground(p) == f


def test_dict_LC_wrt_last():
    R, x, y = ring("x, y", ZZ)
    f = x**2 * y**2 + 3 * x**2 + 5 * x * y + 7
    R_, y_ = ring("y", ZZ)
    expected_lc = y_**2 + 3

    assert dict_LC_wrt_last(dict(f), R.ngens, R.domain) == dict(expected_lc)


def test_dict_deg_wrt_last():
    R, x, y, z = ring("x, y, z", ZZ)
    f1 = x**3 * y**2 * z**5 + x**4 * y * z + x**2 * y**3 * z**8
    assert dict_deg_wrt_last(dict(f1), R.ngens) == (4, 1)

    f2 = x**2 * y * z**3 + x * y**100 * z**5
    assert dict_deg_wrt_last(dict(f2), R.ngens) == (2, 1)


def test_dict_chinese_remainder_reconstruction_multivariate():
    R, x, y = ring("x, y", ZZ)
    p, q = MPZ(3), MPZ(5)

    hp = x**3*y - x**2 - 1
    hq = -x**3*y - 2*x*y**2 + 2

    hpq_dict = dict_chinese_remainder_reconstruction_multivariate(
        dict(hp), dict(hq), p, q, R.domain, R.ngens)
    hpq = R.from_dict(hpq_dict)

    assert hpq.trunc_ground(p) == hp
    assert hpq.trunc_ground(q) == hq

    R, x, y, z = ring("x, y, z", ZZ)
    p, q = MPZ(6), MPZ(5)

    hp = 3*x**4 - y**3*z + z
    hq = -2*x**4 + z

    hpq_dict = dict_chinese_remainder_reconstruction_multivariate(
        dict(hp), dict(hq), p, q, R.domain, R.ngens)
    hpq = R.from_dict(hpq_dict)

    assert hpq.trunc_ground(p) == hp
    assert hpq.trunc_ground(q) == hq
