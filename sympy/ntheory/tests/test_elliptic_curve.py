from __future__ import annotations
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
    # Issue 28529: P + (-P) must return point at infinity on general Weierstrass curves.
    # The curve y^2 + y = x^3 + x^2 has a3=1, so the additive inverse of (x, y)
    # is (x, -y - a3) = (x, -y - 1), not (x, -y). The check y1+y2==0 was wrong.
    e4 = EllipticCurve(0, 0, 0, 1, 1)
    P = e4(0, 0)
    assert (P + (-P)).z == 0
