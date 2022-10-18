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

    e1 = EllipticCurve(-17, 16, modulus=307)
    assert e1(7, 213, 1)

    # Discrete Log
    e1 = EllipticCurve(319, 1007, modulus=2819)
    p1 = e1(1201, 711)
    p2 = e1(506, 435)
    assert p1.discrete_log(p2) == 356
    assert p1*356 == p2
    assert p2.discrete_log(p1) == 421
    assert p2*421 == p1
    e1 = EllipticCurve(321, 231, modulus=509)
    p1 = e1(422, 98)
    p2 = e1(144, 210)
    assert p1.discrete_log(p2) == 215
    assert p1*215 == p2
    assert p2.discrete_log(p1) == 17
    assert p2*17 == p1
