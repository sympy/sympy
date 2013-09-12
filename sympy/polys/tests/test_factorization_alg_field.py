from sympy.polys.factorization_alg_field import efactor
from sympy.polys.domains import AlgebraicField, QQ
from sympy import sqrt, I, ring
from sympy.utilities.pytest import slow

def test_efactor():
    A = AlgebraicField(QQ, sqrt(2))
    R, x, y = ring('x, y', A)

    f = (x**2 + sqrt(2)*y)
    assert efactor(f) == (A.one, [(f, 1)])

    f1 = x + y
    f2 = x**2 + sqrt(2)*y
    f = f1 * f2

    assert efactor(f) == (A.one, [(f1, 1), (f2, 1)])

    f1 = x + 2**10*y
    f2 = x**2 + sqrt(2)*y
    f = f1 * f2

    assert efactor(f) == (A.one, [(f1, 1), (f2, 1)])

    A = AlgebraicField(QQ, sqrt(3))
    R, x, y, z = ring('x, y, z', A)

    f1 = z + 1
    f2 = A([QQ(3, 4)])*x*y**2 + sqrt(3)
    f3 = sqrt(3)*y*x**2 + 2*y + z
    f = f1 * f2**2 * f3

    lc2 = f2.LC
    lc3 = f3.LC

    assert efactor(f) == (lc2**2*lc3, [(f1, 1), (f2.quo_ground(lc2), 2),
        (f3.quo_ground(lc3), 1)])

    A = AlgebraicField(QQ, I)
    R, x, y = ring('x, y', A)

    f1 = x*(y + 1) + 1
    f2 = x*(y + I) + 1
    f3 = x**2*(y - I) + 1
    f = f1*f2*f3

    assert efactor(f) == (A.one, [(f1, 1), (f2, 1), (f3, 1)])

    lc = -A([2])
    f1 = x*(y - 3*I) + lc**(-1)
    f2 = x*(y + I) + 1
    f3 = x*(y + 2) + 1
    f = lc*f1*f2*f3

    assert efactor(f) == (lc, [(f1, 1), (f2, 1), (f3, 1)])

    A = AlgebraicField(QQ, sqrt(8))
    R, x, y = ring('x, y', A)

    f = (x - sqrt(8)/2)*(x + sqrt(8)/2)

    assert efactor(f) == f.factor_list()

    A = AlgebraicField(QQ, sqrt(3))
    R, x, y, z, t = ring('x, y, z, t', A)

    f1 = z*t + 1
    f2 = x*y - sqrt(3)
    f3 = x**2 + 1
    f4 = x**2 + z*t
    f = f1*f2*f3*f4

    assert efactor(f) == (A.one, [(f1, 1), (f2, 1), (f3, 1), (f4, 1)])

    a = QQ(1, 2)*sqrt(2)*(1 + I)
    A = AlgebraicField(QQ, a)
    R, x, y = ring('x, y', A)

    f = x**4 + y**4

    f1 = x - a*y
    f2 = x - a**3*y
    f3 = x + a*y
    f4 = x + a**3*y

    assert efactor(f) == (A.one, [(f1, 1), (f2, 1), (f3, 1), (f4, 1)])


@slow
def test_efactor_wang():
    a = QQ(1, 4)*(-1 + sqrt(5)) - I*sqrt(QQ(1, 8)*(sqrt(5) + 5))
    A = AlgebraicField(QQ, a)
    R, x, y, z = ring('x, y, z', A)

    f = x**8 + 2*x**7 - (y + z**2 + 8)*x**6 + 2*(-2*y + 3*z**2 - 20)*x**5 + \
        (y**2 + 2*(z**2 - 24)*y + z**4 + 32*z**2 + 256)*x**4 + \
        (-4*y**2 + 2*(z**2 + 16)*y - 4*z**4 + 32*z**2 + 960)*x**3 + \
        (-y**3 + (-3*z**2 + 28)*y**2 + 2*(z**4 - 2*z**2 + 192)*y - z**6 -
        32*z**4 + 144*z**2 - 1152)*x**2 + (2*y**3 + 4*(-z**2 + 18)*y**2 +
        6*(z**4 + 4*z**2 - 96)*y + 2*z**6 - 48*z**4 - 576*z**2 + 3456)*x + \
        y**4 - (z**2 + 12)*y**3 + (z**4 + 24*z**2 + 144)*y**2 - \
        (z**6 - 24*z**4 + 432*z**2 + 1728)*y + z**8 - 12*z**6 + 144*z**4 - \
        1728*z**2 + 20736

    f1 = x**2 - 2*a*x - (a**3 + a**2 + a + 1)*z**2 + a**2*y + 12*a**3
    f2 = x**2 - 2*a**2*x + a**3*z**2 - (a**3 + a**2 + a + 1)*y + 12*a
    f3 = x**2 - 2*a**3*x + a**2*z**2 + a*y - 12*(a**3 + a**2 + a + 1)
    f4 = x**2 + 2*(a**3 + a**2 + a + 1)*x + a*z**2 + a**3*y + 12*a**2

    assert efactor(f) == (A.one, [(f1, 1), (f2, 1), (f3, 1), (f4, 1)])
