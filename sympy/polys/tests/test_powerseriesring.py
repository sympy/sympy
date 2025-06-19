from sympy import QQ, ZZ
from sympy.polys.powerseriesring import power_series_ring


def test_arithematic_zz():
    R = power_series_ring(ZZ)
    x = R.gen
    assert R.pretty(x) == 'x'

    p1 = R.multiply(R.add(x, R.one), x)
    s1 = R.multiply(p1, R.multiply(p1, p1))
    s2 = R.subtract(s1, p1)

    assert R.pretty(p1) == 'x**2 + x'
    assert R.pretty(s1) == 'x**3 + 3*x**4 + 3*x**5 + O(x**6)'
    assert R.pretty(s2) == '-x - x**2 + x**3 + 3*x**4 + 3*x**5 + O(x**6)'

    p2 = R.multiply_ground(R.one, 5)
    p3 = R.add(p2, R.multiply_ground(x, -3))
    p4 = R.multiply(p3, p3)

    assert R.pretty(p2) == '5'
    assert R.pretty(p3) == '-3*x + 5'
    assert R.pretty(p4) == '9*x**2 - 30*x + 25'

    p5 = R.subtract(R.zero, x)
    p6 = R.multiply(p5, R.multiply(p5, p5))

    assert R.pretty(p5) == '-x'
    assert R.pretty(p6) == '-x**3'

    p7 = R.add(R.one, R.multiply_ground(R.multiply(x, x), 2))
    p8 = R.subtract(p7, R.multiply_ground(R.multiply(R.multiply(x, x), x), -4))
    p9 = R.multiply_ground(R.multiply(R.multiply(x, x), R.multiply(x, x)), 7)
    p10 = R.add(p8, p9)

    assert R.pretty(p7) == '2*x**2 + 1'
    assert R.pretty(p8) == '4*x**3 + 2*x**2 + 1'
    assert R.pretty(p9) == '7*x**4'
    assert R.pretty(p10) == '7*x**4 + 4*x**3 + 2*x**2 + 1'

    p11 = R.multiply_ground(R.multiply(x, x), -1)
    p12 = R.add(p11, R.multiply_ground(R.multiply(R.multiply(x, x), x), 6))
    p13 = R.multiply(p12, x)

    assert R.pretty(p11) == '-x**2'
    assert R.pretty(p12) == '6*x**3 - x**2'
    assert R.pretty(p13) == '6*x**4 - x**3'

def test_arithematic_qq():
    R = power_series_ring(QQ)
    x = R.gen
    assert R.pretty(x) == 'x'

    p1 = R.multiply_ground(R.one, QQ(1, 2))
    p2 = R.add(p1, R.multiply_ground(x, QQ(3, 4)))
    p3 = R.multiply(p2, p2)

    assert R.pretty(p1) == '1/2'
    assert R.pretty(p2) == '3/4*x + 1/2'
    assert R.pretty(p3) == '9/16*x**2 + 3/4*x + 1/4'

    p4 = R.subtract(R.one, R.multiply_ground(x, QQ(1, 3)))
    p5 = R.multiply(p4, R.multiply(p4, x))

    assert R.pretty(p4) == '-1/3*x + 1'
    assert R.pretty(p5) == '1/9*x**3 - 2/3*x**2 + x'

    p6 = R.add(R.multiply_ground(R.one, QQ(-2, 5)), R.multiply_ground(R.multiply(x, x), QQ(7, 8)))
    p7 = R.subtract(p6, R.multiply_ground(R.multiply(R.multiply(x, x), x), QQ(1, 6)))

    assert R.pretty(p6) == '7/8*x**2 - 2/5'
    assert R.pretty(p7) == '-1/6*x**3 + 7/8*x**2 - 2/5'

    p8 = R.multiply_ground(x, QQ(-1, 4))
    p9 = R.multiply(p8, R.multiply(p8, p8))

    assert R.pretty(p8) == '-1/4*x'
    assert R.pretty(p9) == '-1/64*x**3'

    p10 = R.add(R.multiply_ground(R.one, QQ(3, 7)), R.multiply_ground(R.multiply(x, x), QQ(-5, 9)))
    p11 = R.multiply(p10, R.multiply(x, x))
    p12 = R.add(p11, R.multiply_ground(R.multiply(R.multiply(x, x), R.multiply(x, x)), QQ(2, 11)))

    assert R.pretty(p10) == '-5/9*x**2 + 3/7'
    assert R.pretty(p11) == '-5/9*x**4 + 3/7*x**2'
    assert R.pretty(p12) == '-37/99*x**4 + 3/7*x**2'

    p13 = R.multiply_ground(R.multiply(R.multiply(x, x), x), QQ(-4, 13))
    p14 = R.subtract(R.multiply_ground(R.one, QQ(1, 2)), p13)
    p15 = R.multiply(p14, p14)

    assert R.pretty(p13) == '-4/13*x**3'
    assert R.pretty(p14) == '4/13*x**3 + 1/2'
    assert R.pretty(p15) == '1/4 + 4/13*x**3 + O(x**6)'
