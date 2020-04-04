from sympy.ntheory.elliptic_curve import EllipticCurve
from sympy.polys import Poly

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

    e2 = EllipticCurve(0, 1, modulus=3)
    X = e2.X
    Y = e2.Y
    assert e2.div_poly(6) == Poly(0, X, Y, modulus=3)
    assert e2.div_poly(7) == Poly(X**18*Y**4 - X**9*Y**4 + Y**4, X, Y, modulus=3)
