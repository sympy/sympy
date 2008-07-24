from sympy import Symbol, zeta, nan, Rational, Real, pi, dirichlet_eta, log, zoo

x = Symbol('x')

def test_zeta():

    assert zeta(nan) == nan
    assert zeta(x, nan) == nan

    assert zeta(0) == Rational(-1,2)
    assert zeta(0, x) == Rational(1,2) - x

    assert zeta(1) == zoo
    assert zeta(1, 2) == zoo
    assert zeta(1, -7) == zoo
    assert zeta(1, x) == zoo

    assert zeta(2, 0) == pi**2/6
    assert zeta(2, 1) == pi**2/6

    assert zeta(2) == pi**2/6
    assert zeta(4) == pi**4/90
    assert zeta(6) == pi**6/945

    assert zeta(2, 2) == pi**2/6 - 1
    assert zeta(4, 3) == pi**4/90 - Rational(17, 16)
    assert zeta(6, 4) == pi**6/945 - Rational(47449, 46656)

    assert zeta(2, -2) == pi**2/6 + Rational(5, 4)
    assert zeta(4, -3) == pi**4/90 + Rational(1393, 1296)
    assert zeta(6, -4) == pi**6/945 + Rational(3037465, 2985984)

    assert zeta(-1) == -Rational(1, 12)
    assert zeta(-2) == 0
    assert zeta(-3) == Rational(1, 120)
    assert zeta(-4) == 0
    assert zeta(-5) == -Rational(1, 252)

    assert zeta(-1, 3) == -Rational(37, 12)
    assert zeta(-1, 7) == -Rational(253, 12)
    assert zeta(-1, -4) == Rational(119, 12)
    assert zeta(-1, -9) == Rational(539, 12)

    assert zeta(-4, 3) == -17
    assert zeta(-4, -8) == 8772

    assert zeta(0, 0) == -Rational(1, 2)

    assert zeta(0, 1) == -Rational(1, 2)
    assert zeta(0, -1) == Rational(1, 2)

    assert zeta(0, 2) == -Rational(3, 2)
    assert zeta(0, -2) == Rational(3, 2)

    assert zeta(3).evalf(20).epsilon_eq(Real("1.2020569031595942854",20), 1e-19)

def test_dirichlet_eta():

    assert dirichlet_eta(0) == Rational(1,2)
    assert dirichlet_eta(-1) == Rational(1,4)
    assert dirichlet_eta(1) == log(2)
    assert dirichlet_eta(2) == pi**2/12
    assert dirichlet_eta(4) == pi**4*Rational(7,720)
