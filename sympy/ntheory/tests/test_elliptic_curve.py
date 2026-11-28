from __future__ import annotations
from sympy import Rational
from sympy.ntheory.elliptic_curve import EllipticCurve, EllipticCurvePoint


def test_elliptic_curve_discriminant():
    curve = EllipticCurve(Rational(1, 8), 0)
    assert curve.discriminant == Rational(-1, 8)
    assert not curve.is_singular

    assert EllipticCurve(0, 0).is_singular

    assert type(EllipticCurve(1, 0).discriminant) is int
    finite_field_curve = EllipticCurve(1, 0, modulus=5)
    assert type(finite_field_curve.discriminant) is int
    assert finite_field_curve.discriminant == 1
    assert not finite_field_curve.is_singular
    assert EllipticCurve(0, 0, modulus=5).is_singular


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
