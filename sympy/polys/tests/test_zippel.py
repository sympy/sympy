from __future__ import annotations
from sympy.polys.rings import ring
from sympy.polys.domains import ZZ
from sympy.polys.zippel import (_LC, _chinese_remainder_reconstruction_multivariate, _deg,
    _gf_gcd, _trivial_gcd, _primitive, lag_basis, vandermonde_interp)
from sympy.matrices import Matrix


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
