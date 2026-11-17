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


def test_elliptic_curve_point_order():
    # Issue 31160: EllipticCurvePoint.order() correctness
    # Point at infinity has order 1 and acts as additive identity
    e_q = EllipticCurve(-3, 3)
    O = EllipticCurvePoint.point_at_infinity(e_q)
    assert O.order() == 1

    # Short Weierstrass over QQ: horizontal tangent point has (2*P).y == -P.y
    # with (2*P).x != P.x; infinite order must not be misclassified as order 3
    p_inf = e_q(1, 1)
    assert (p_inf + O).x == p_inf.x and (p_inf + O).y == p_inf.y
    assert (O + p_inf).x == p_inf.x and (O + p_inf).y == p_inf.y
    assert p_inf.order() == oo
    assert EllipticCurve(-1, 9)(0, 3).order() == oo
    tors = e_q.torsion_points()
    assert len(tors) == 1 and tors[0].z == 0

    # Order-2 point on short Weierstrass model
    p_short2 = EllipticCurve(0, 1)(-1, 0)
    assert p_short2.order() == 2
    assert (2 * p_short2).z == 0 and (1 * p_short2).z != 0

    # Order-2 points with y != 0 on general Weierstrass models
    for pt in (EllipticCurve(1, -1, a3=2)(0, -1), EllipticCurve(1, -4, a3=4)(0, -2)):
        assert pt.order() == 2
        assert (2 * pt).z == 0 and (1 * pt).z != 0

    # Order-3 points on short and general Weierstrass models
    for pt in (EllipticCurve(0, 1)(0, 1), EllipticCurve(0, 0, a3=1)(0, 0)):
        assert pt.order() == 3
        assert (3 * pt).z == 0 and (2 * pt).z != 0

    # General Weierstrass over QQ: y == 0 does not imply order 2
    p_gen5 = EllipticCurve(0, 0, a1=0, a2=-1, a3=1)(0, 0)
    assert p_gen5.order() == 5
    assert (5 * p_gen5).z == 0 and (4 * p_gen5).z != 0

    # Higher torsion orders over QQ (including rational coefficients and order 12)
    for pt, expected in (
        (EllipticCurve(-3, 3, a1=1, a2=-1, a3=1)(1, 0), 7),
        (EllipticCurve(0, 0, a1=Rational(-1, 2), a2=-3, a3=-3)(0, 0), 8),
        (EllipticCurve(0, 0, a1=43, a2=-210, a3=-210)(0, 0), 12),
    ):
        assert pt.order() == expected
        assert (expected * pt).z == 0
        assert all((k * pt).z != 0 for k in range(1, expected))

    # Rational-coefficient curve with fractional torsion coordinates and z != 1
    e_rat = EllipticCurve(0, Rational(1, 64))
    for pt in (
        e_rat(Rational(1, 2), Rational(3, 8)),
        e_rat(Rational(3, 2), Rational(9, 8), 3),
    ):
        assert pt.order() == 6
        assert (6 * pt).z == 0 and (5 * pt).z != 0

    # Short Weierstrass over FF(5): (2*P).y == -P.y with (2*P).x != P.x
    p_fp = EllipticCurve(2, 1, modulus=5)(1, 2)
    assert p_fp.order() == 7
    assert (7 * p_fp).z == 0 and (6 * p_fp).z != 0

    # Characteristic 2 finite field
    e_f2 = EllipticCurve(0, 1, a1=1, a3=1, modulus=2)
    for x, y in e_f2.points():
        pt = e_f2(x, y)
        n = pt.order()
        assert (n * pt).z == 0
        assert all((k * pt).z != 0 for k in range(1, n))
