from __future__ import annotations
from sympy.core.numbers import Rational, oo
from sympy.ntheory.elliptic_curve import EllipticCurve, EllipticCurvePoint


def test_elliptic_curve():
    # Point addition and multiplication
    e3 = EllipticCurve(-1, 9)
    p = e3(0, 3)
    q = e3(-1, 3)
    r = p + q
    assert r.x == 1 and r.y == -3
    r = 2*p + q
    assert r.x == 35 and r.y == 207
    r = -p + q
    assert r.x == 37 and r.y == 225
    # Verify result in http://www.lmfdb.org/EllipticCurve/Q
    # Discriminant
    assert EllipticCurve(-1, 9).discriminant == -34928
    assert EllipticCurve(-2731, -55146, 1, 0, 1).discriminant == 25088
    # Torsion points
    assert len(EllipticCurve(0, 1).torsion_points()) == 6
    # Issue 28546: -O should return canonical infinity point
    O = EllipticCurvePoint.point_at_infinity(e3)
    assert (-O).x == O.x and (-O).y == O.y and (-O).z == O.z
    # Issue 28529
    curv = EllipticCurve(0, 0, 0, 1, 1, 0)
    assert (curv(0,0,1) + -curv(0,0,1)).z == 0


def test_elliptic_curve_point_order_infinity():
    curve = EllipticCurve(-3, 3)
    O = EllipticCurvePoint.point_at_infinity(curve)
    pt = curve(1, 1)
    assert O.order() == 1
    assert (pt + O).x == pt.x and (pt + O).y == pt.y and (pt + O).z == pt.z
    assert (O + pt).x == pt.x and (O + pt).y == pt.y and (O + pt).z == pt.z


def test_elliptic_curve_point_order_2():
    # Short Weierstrass order-2 point (y == 0)
    pt_short = EllipticCurve(0, 1)(-1, 0)
    assert pt_short.order() == 2
    assert (2 * pt_short).z == 0 and (1 * pt_short).z != 0

    # General Weierstrass order-2 points with y != 0 (false negatives previously)
    for pt in (
        EllipticCurve(1, -1, a3=2)(0, -1),
        EllipticCurve(1, -4, a3=4)(0, -2),
    ):
        assert pt.order() == 2
        assert (2 * pt).z == 0 and (1 * pt).z != 0

    # General Weierstrass y == 0 point of order 5 (false positive previously)
    pt_gen5 = EllipticCurve(0, 0, a1=0, a2=-1, a3=1)(0, 0)
    assert pt_gen5.order() == 5
    assert (5 * pt_gen5).z == 0 and (4 * pt_gen5).z != 0


def test_elliptic_curve_point_order_3():
    # Genuine order-3 points on short and general Weierstrass models
    for pt in (
        EllipticCurve(0, 1)(0, 1),
        EllipticCurve(0, 0, a3=1)(0, 0),
    ):
        assert pt.order() == 3
        assert (3 * pt).z == 0 and (2 * pt).z != 0

    # Horizontal-tangent point: (2*P).y == -P.y with (2*P).x != P.x (not order 3)
    pt_inf = EllipticCurve(-3, 3)(1, 1)
    assert (2 * pt_inf).y == -pt_inf.y
    assert (3 * pt_inf).z != 0
    assert pt_inf.order() == oo


def test_elliptic_curve_point_order_rational():
    # Fractional-coordinate order-6 point (z == 1 and projective z == 3)
    curve_rat = EllipticCurve(0, Rational(1, 64))
    for pt in (
        curve_rat(Rational(1, 2), Rational(3, 8)),
        curve_rat(Rational(3, 2), Rational(9, 8), 3),
    ):
        assert pt.order() == 6
        assert (6 * pt).z == 0 and (5 * pt).z != 0

    # Higher rational torsion orders including rational coefficients and order 12
    for pt, expected in (
        (EllipticCurve(-3, 3, a1=1, a2=-1, a3=1)(1, 0), 7),
        (EllipticCurve(0, 0, a1=Rational(-1, 2), a2=-3, a3=-3)(0, 0), 8),
        (EllipticCurve(0, 0, a1=43, a2=-210, a3=-210)(0, 0), 12),
    ):
        assert pt.order() == expected
        assert (expected * pt).z == 0
        assert all((k * pt).z != 0 for k in range(1, expected))


def test_elliptic_curve_point_order_infinite():
    assert EllipticCurve(-3, 3)(1, 1).order() == oo
    assert EllipticCurve(-1, 9)(0, 3).order() == oo


def test_elliptic_curve_point_order_finite_field():
    # FF(5) order-7 point where (2*P).y == -P.y with (2*P).x != P.x
    pt_fp = EllipticCurve(2, 1, modulus=5)(1, 2)
    assert pt_fp.order() == 7
    assert (7 * pt_fp).z == 0 and (6 * pt_fp).z != 0

    # Characteristic-2 finite field: verify minimal order for every point
    curve_f2 = EllipticCurve(0, 1, a1=1, a3=1, modulus=2)
    for x, y in curve_f2.points():
        pt = curve_f2(x, y)
        n = pt.order()
        assert (n * pt).z == 0
        assert all((k * pt).z != 0 for k in range(1, n))


def test_elliptic_curve_torsion_points_order_regression():
    # Issue 31160: infinite-order points (1, 1) and (1, -1) must not appear
    tors = EllipticCurve(-3, 3).torsion_points()
    assert len(tors) == 1 and tors[0].z == 0
