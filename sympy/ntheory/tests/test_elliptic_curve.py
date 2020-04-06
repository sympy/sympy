from sympy.ntheory.elliptic_curve import EllipticCurve

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

    assert EllipticCurve(46, 74, modulus=101).schoof() == 97
    assert EllipticCurve(101 ,210 , modulus=317).schoof() == 308
    assert EllipticCurve(132 ,132 , modulus=401).schoof() == 383
