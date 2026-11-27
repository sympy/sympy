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


def test_point_order():
    # Issue 31160
    e = EllipticCurve(-3, 3)
    assert e(1, 1).order() == oo
    assert len(e.torsion_points()) == 1
    assert EllipticCurve(2, 1, modulus=5)(1, 2).order() == 7
    assert EllipticCurve(0, 0, 0, -1, 1)(0, 0).order() == 5
    assert EllipticCurve(-1, 0, 0, 0, 1)(0, 0).order() == oo
    assert EllipticCurve(1, -1, a3=2)(0, -1).order() == 2
    assert EllipticCurve(1, -4, a3=4)(0, -2).order() == 2
    assert EllipticCurve(-1, -1, -3, 0, 1)(Rational(-1, 4), Rational(-7, 8)).order() == 2
    assert EllipticCurve(0, Rational(1, 64))(Rational(1, 2), Rational(3, 8)).order() == 6
    # Tate normal form y**2 + (1 - c)*x*y - b*y = x**3 - b*x**2
    for n, b, c in [(4, 2, 0), (5, 2, 2), (6, 6, 2), (7, 4, 2),
                    (8, 3, Rational(3, 2)), (9, 12, 4), (10, 24, 6),
                    (12, 210, -42)]:
        assert EllipticCurve(0, 0, 1 - c, -b, -b)(0, 0).order() == n
    e = EllipticCurve(4, 6, 2, 1, 6, modulus=7)
    for x, y in e.points():
        p = e(x, y)
        n = p.order()
        assert (n*p).z == 0
        assert all((m*p).z != 0 for m in range(1, n))
