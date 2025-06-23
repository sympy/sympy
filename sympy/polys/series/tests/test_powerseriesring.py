from sympy.polys.domains import QQ, ZZ
from sympy.polys.series.seriesring import power_series_ring


def test_edge_cases():
    R = power_series_ring(ZZ, 0)
    x = R.gen
    assert R.pretty(x) == 'O(x**0)'
    assert R.pretty(R.zero) == 'O(x**0)'
    assert R.pretty(R.one) == 'O(x**0)'

    R = power_series_ring(QQ, 1)
    x = R.gen
    assert R.pretty(x) == '0 + O(x**1)'
    assert R.pretty(R.zero) == '0'
    assert R.pretty(R.one) == '1'


def test_equality():

    # Ring equality
    R1 = power_series_ring(ZZ)
    R2 = power_series_ring(QQ, 5)

    assert R1 == power_series_ring(ZZ, prec=6)
    assert R1 != R2

    # Series equality
    x = R1.gen
    x2 = R2.gen

    assert x == x2

    s1 = R1.pow_int(R1.multiply(R1.add(x, R1.one), x), 6)
    s2 = R2.pow_int(R2.multiply(x2, x2), 6)

    assert R1.equal(s1, s1) is None
    assert R1.equal(s1, s1) is None
    assert R1.equal(s1, s2) is False


def test_arithematic_zz():
    R = power_series_ring(ZZ)
    x = R.gen
    assert R.pretty(x) == 'x'

    p1 = R.multiply(R.add(x, R.one), x)
    s1 = R.multiply(p1, R.multiply(p1, p1))
    s2 = R.subtract(s1, p1)

    assert R.pretty(p1) == 'x + x**2'
    assert R.pretty(s1) == 'x**3 + 3*x**4 + 3*x**5 + O(x**6)'
    assert R.pretty(s2) == '-x - x**2 + x**3 + 3*x**4 + 3*x**5 + O(x**6)'

    p2 = R.multiply_ground(R.one, 5)
    p3 = R.add(p2, R.multiply_ground(x, -3))
    p4 = R.multiply(p3, p3)

    assert R.pretty(p2) == '5'
    assert R.pretty(p3) == '5 - 3*x'
    assert R.pretty(p4) == '25 - 30*x + 9*x**2'

    p5 = R.subtract(R.zero, x)
    p6 = R.multiply(p5, R.multiply(p5, p5))

    assert R.pretty(p5) == '-x'
    assert R.pretty(p6) == '-x**3'

    p7 = R.add(R.one, R.multiply_ground(R.multiply(x, x), 2))
    p8 = R.subtract(p7, R.multiply_ground(R.multiply(R.multiply(x, x), x),
                                          -4))
    p9 = R.multiply_ground(R.multiply(R.multiply(x, x),
                                      R.multiply(x, x)), 7)
    p10 = R.add(p8, p9)

    assert R.pretty(p7) == '1 + 2*x**2'
    assert R.pretty(p8) == '1 + 2*x**2 + 4*x**3'
    assert R.pretty(p9) == '7*x**4'
    assert R.pretty(p10) == '1 + 2*x**2 + 4*x**3 + 7*x**4'

    p11 = R.multiply_ground(R.multiply(x, x), -1)
    p12 = R.add(p11, R.multiply_ground(R.multiply(R.multiply(x, x), x), 6))
    p13 = R.multiply(p12, x)

    assert R.pretty(p11) == '-x**2'
    assert R.pretty(p12) == '-x**2 + 6*x**3'
    assert R.pretty(p13) == '-x**3 + 6*x**4'

    sq1 = R.square(x)
    sq2 = R.square(R.add(x, R.one))
    sq3 = R.square(R.multiply_ground(x, 3))
    sq4 = R.square(R.add(R.multiply_ground(x, 2),
                         R.multiply_ground(R.one, 3)))
    sq5 = R.square(R.subtract(R.multiply_ground(x, 4),
                              R.multiply_ground(R.one, 1)))
    sq6 = R.square(R.add(R.multiply_ground(R.multiply(x, x), 2),
                         R.multiply_ground(x, -5)))

    assert R.pretty(sq1) == 'x**2'
    assert R.pretty(sq2) == '1 + 2*x + x**2'
    assert R.pretty(sq3) == '9*x**2'
    assert R.pretty(sq4) == '9 + 12*x + 4*x**2'
    assert R.pretty(sq5) == '1 - 8*x + 16*x**2'
    assert R.pretty(sq6) == '25*x**2 - 20*x**3 + 4*x**4'

    pow1 = R.pow_int(x, 3)
    pow2 = R.pow_int(R.add(x, R.one), 2)
    pow3 = R.pow_int(R.multiply_ground(x, 2), 4)
    pow4 = R.pow_int(R.add(R.multiply_ground(x, 2), R.one), 7)

    assert R.pretty(pow1) == 'x**3'
    assert R.pretty(pow2) == '1 + 2*x + x**2'
    assert R.pretty(pow3) == '16*x**4'
    assert R.pretty(pow4) == ('1 + 14*x + 84*x**2 + 280*x**3 + 560*x**4 + '
                              '672*x**5 + O(x**6)')


def test_arithematic_qq():
    R = power_series_ring(QQ)
    x = R.gen
    assert R.pretty(x) == 'x'

    p1 = R.multiply_ground(R.one, QQ(1, 2))
    p2 = R.add(p1, R.multiply_ground(x, QQ(3, 4)))
    p3 = R.multiply(p2, p2)

    assert R.pretty(p1) == '1/2'
    assert R.pretty(p2) == '1/2 + 3/4*x'
    assert R.pretty(p3) == '1/4 + 3/4*x + 9/16*x**2'

    p4 = R.subtract(R.one, R.multiply_ground(x, QQ(1, 3)))
    p5 = R.multiply(p4, R.multiply(p4, x))

    assert R.pretty(p4) == '1 - 1/3*x'
    assert R.pretty(p5) == 'x - 2/3*x**2 + 1/9*x**3'

    p6 = R.add(R.multiply_ground(R.one, QQ(-2, 5)),
                R.multiply_ground(R.multiply(x, x), QQ(7, 8)))
    p7 = R.subtract(p6, R.multiply_ground(R.multiply(R.multiply(x, x), x),
                                          QQ(1, 6)))

    assert R.pretty(p6) == '-2/5 + 7/8*x**2'
    assert R.pretty(p7) == '-2/5 + 7/8*x**2 - 1/6*x**3'

    p8 = R.multiply_ground(x, QQ(-1, 4))
    p9 = R.multiply(p8, R.multiply(p8, p8))

    assert R.pretty(p8) == '-1/4*x'
    assert R.pretty(p9) == '-1/64*x**3'

    p10 = R.add(R.multiply_ground(R.one, QQ(3, 7)),
                 R.multiply_ground(R.multiply(x, x), QQ(-5, 9)))
    p11 = R.multiply(p10, R.multiply(x, x))
    p12 = R.add(p11, R.multiply_ground(R.multiply(R.multiply(x, x),
                                                  R.multiply(x, x)),
                                       QQ(2, 11)))

    assert R.pretty(p10) == '3/7 - 5/9*x**2'
    assert R.pretty(p11) == '3/7*x**2 - 5/9*x**4'
    assert R.pretty(p12) == '3/7*x**2 - 37/99*x**4'

    p13 = R.multiply_ground(R.multiply(R.multiply(x, x), x), QQ(-4, 13))
    p14 = R.subtract(R.multiply_ground(R.one, QQ(1, 2)), p13)
    p15 = R.multiply(p14, p14)

    assert R.pretty(p13) == '-4/13*x**3'
    assert R.pretty(p14) == '1/2 + 4/13*x**3'
    assert R.pretty(p15) == '1/4 + 4/13*x**3 + O(x**6)'

    sq1 = R.square(x)
    sq2 = R.square(R.add(x, R.multiply_ground(R.one, QQ(1, 2))))
    sq3 = R.square(R.multiply_ground(x, QQ(2, 3)))
    sq4 = R.square(R.add(R.multiply_ground(x, QQ(3, 4)),
                         R.multiply_ground(R.one, QQ(2, 5))))
    sq5 = R.square(R.subtract(R.multiply_ground(x, QQ(5, 6)),
                              R.multiply_ground(R.one, QQ(1, 7))))
    sq6 = R.square(R.add(R.multiply_ground(R.multiply(x, x), QQ(1, 3)),
                         R.multiply_ground(x, QQ(-2, 5))))

    assert R.pretty(sq1) == 'x**2'
    assert R.pretty(sq2) == '1/4 + x + x**2'
    assert R.pretty(sq3) == '4/9*x**2'
    assert R.pretty(sq4) == '4/25 + 3/5*x + 9/16*x**2'
    assert R.pretty(sq5) == '1/49 - 5/21*x + 25/36*x**2'
    assert R.pretty(sq6) == '4/25*x**2 - 4/15*x**3 + 1/9*x**4'

    pow1 = R.pow_int(x, 11)
    pow2 = R.pow_int(R.add(x, R.multiply_ground(R.one, QQ(1, 3))), 3)
    pow3 = R.pow_int(R.multiply_ground(x, QQ(3, 2)), 4)
    pow4 = R.pow_int(R.add(R.multiply_ground(x, QQ(1, 2)),
                           R.multiply_ground(R.one, QQ(2, 3))), 4)

    assert R.pretty(pow1) == '0 + O(x**6)'
    assert R.pretty(pow2) == '1/27 + 1/3*x + x**2 + x**3'
    assert R.pretty(pow3) == '81/16*x**4'
    assert R.pretty(pow4) == ('16/81 + 16/27*x + 2/3*x**2 + 1/3*x**3 + '
                              '1/16*x**4')
