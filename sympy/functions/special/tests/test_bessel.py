from sympy import jn, yn, symbols, sin, cos, pi, S, jn_zeros, besselj, \
    bessely, besseli, besselk, hankel1, hankel2, expand_func, \
    latex, sqrt, sinh, cosh
from sympy.functions.special.bessel import fn
from sympy.utilities.pytest import raises, skip
from sympy.utilities.randtest import \
    random_complex_number as randcplx, \
    test_numerically as tn, \
    test_derivative_numerically as td, \
    _randint

from sympy.abc import z, n, k, x

randint = _randint()

def test_bessel_rand():
    for f in [besselj, bessely, besseli, besselk, hankel1, hankel2, jn, yn]:
        assert td(f(randcplx(), z), z)


def test_diff():
    assert besselj(n, z).diff(z) == besselj(n - 1, z)/2 - besselj(n + 1, z)/2
    assert bessely(n, z).diff(z) == bessely(n - 1, z)/2 - bessely(n + 1, z)/2
    assert besseli(n, z).diff(z) == besseli(n - 1, z)/2 + besseli(n + 1, z)/2
    assert besselk(n, z).diff(z) == -besselk(n - 1, z)/2 - besselk(n + 1, z)/2
    assert hankel1(n, z).diff(z) == hankel1(n - 1, z)/2 - hankel1(n + 1, z)/2
    assert hankel2(n, z).diff(z) == hankel2(n - 1, z)/2 - hankel2(n + 1, z)/2


def test_rewrite():
    from sympy import polar_lift, exp, I

    assert besselj(n, z).rewrite(jn) == sqrt(2*z/pi)*jn(n - S(1)/2, z)
    assert bessely(n, z).rewrite(yn) == sqrt(2*z/pi)*yn(n - S(1)/2, z)
    assert besseli(n, z).rewrite(besselj) == \
        exp(-I*n*pi/2)*besselj(n, polar_lift(I)*z)
    assert besselj(n, z).rewrite(besseli) == \
        exp(I*n*pi/2)*besseli(n, polar_lift(-I)*z)

    nu = randcplx()

    assert tn(besselj(nu, z), besselj(nu, z).rewrite(besseli), z)
    assert tn(besselj(nu, z), besselj(nu, z).rewrite(bessely), z)

    assert tn(besseli(nu, z), besseli(nu, z).rewrite(besselj), z)
    assert tn(besseli(nu, z), besseli(nu, z).rewrite(bessely), z)

    assert tn(bessely(nu, z), bessely(nu, z).rewrite(besselj), z)
    assert tn(bessely(nu, z), bessely(nu, z).rewrite(besseli), z)

    assert tn(besselk(nu, z), besselk(nu, z).rewrite(besselj), z)
    assert tn(besselk(nu, z), besselk(nu, z).rewrite(besseli), z)
    assert tn(besselk(nu, z), besselk(nu, z).rewrite(bessely), z)


def test_expand():
    from sympy import besselsimp, Symbol, exp, exp_polar, I

    assert expand_func(besselj(S(1)/2, z).rewrite(jn)) == \
        sqrt(2)*sin(z)/(sqrt(pi)*sqrt(z))
    assert expand_func(bessely(S(1)/2, z).rewrite(yn)) == \
        -sqrt(2)*cos(z)/(sqrt(pi)*sqrt(z))

    # XXX: teach sin/cos to work around arguments like
    # x*exp_polar(I*pi*n/2).  Then change besselsimp -> expand_func
    assert besselsimp(besselj(S(1)/2, z)) ==  sqrt(2)*sin(z)/(sqrt(pi)*sqrt(z))
    assert besselsimp(besselj(S(-1)/2, z)) == sqrt(2)*cos(z)/(sqrt(pi)*sqrt(z))
    assert besselsimp(besselj(S(5)/2, z)) == \
        -sqrt(2)*(z**2*sin(z) + 3*z*cos(z) - 3*sin(z))/(sqrt(pi)*z**(S(5)/2))
    assert besselsimp(besselj(-S(5)/2, z)) == \
        -sqrt(2)*(z**2*cos(z) - 3*z*sin(z) - 3*cos(z))/(sqrt(pi)*z**(S(5)/2))

    assert besselsimp(bessely(S(1)/2, z)) == \
        -(sqrt(2)*cos(z))/(sqrt(pi)*sqrt(z))
    assert besselsimp(bessely(S(-1)/2, z)) == sqrt(2)*sin(z)/(sqrt(pi)*sqrt(z))
    assert besselsimp(bessely(S(5)/2, z)) == \
        sqrt(2)*(z**2*cos(z) - 3*z*sin(z) - 3*cos(z))/(sqrt(pi)*z**(S(5)/2))
    assert besselsimp(bessely(S(-5)/2, z)) == \
        -sqrt(2)*(z**2*sin(z) + 3*z*cos(z) - 3*sin(z))/(sqrt(pi)*z**(S(5)/2))

    assert besselsimp(besseli(S(1)/2, z)) == sqrt(2)*sinh(z)/(sqrt(pi)*sqrt(z))
    assert besselsimp(besseli(S(-1)/2, z)) == \
        sqrt(2)*cosh(z)/(sqrt(pi)*sqrt(z))
    assert besselsimp(besseli(S(5)/2, z)) == \
        sqrt(2)*(z**2*sinh(z) - 3*z*cosh(z) + 3*sinh(z))/(sqrt(pi)*z**(S(5)/2))
    assert besselsimp(besseli(S(-5)/2, z)) == \
        sqrt(2)*(z**2*cosh(z) - 3*z*sinh(z) + 3*cosh(z))/(sqrt(pi)*z**(S(5)/2))

    assert besselsimp(besselk(S(1)/2, z)) == \
        besselsimp(besselk(S(-1)/2, z)) == sqrt(pi)*exp(-z)/(sqrt(2)*sqrt(z))
    assert besselsimp(besselk(S(5)/2, z)) == \
        besselsimp(besselk(S(-5)/2, z)) == \
        sqrt(2)*sqrt(pi)*(z**2 + 3*z + 3)*exp(-z)/(2*z**(S(5)/2))

    def check(eq, ans):
        return tn(eq, ans) and eq == ans

    rn = randcplx(a=1, b=0, d=0, c=2)

    for besselx in [besselj, bessely, besseli, besselk]:
        ri = S(2*randint(-11, 10) + 1) / 2  # half integer in [-21/2, 21/2]
        assert tn(besselsimp(besselx(ri, z)), besselx(ri, z))

    assert check(expand_func(besseli(rn, x)), \
        besseli(rn - 2, x) - 2*(rn - 1)*besseli(rn - 1, x)/x)
    assert check(expand_func(besseli(-rn, x)), \
        besseli(-rn + 2, x) + 2*(-rn + 1)*besseli(-rn + 1, x)/x)

    assert check(expand_func(besselj(rn, x)), \
        -besselj(rn - 2, x) + 2*(rn - 1)*besselj(rn - 1, x)/x)
    assert check(expand_func(besselj(-rn, x)), \
        -besselj(-rn + 2, x) + 2*(-rn + 1)*besselj(-rn + 1, x)/x)

    assert check(expand_func(besselk(rn, x)), \
        besselk(rn - 2, x) + 2*(rn - 1)*besselk(rn - 1, x)/x)
    assert check(expand_func(besselk(-rn, x)), \
        besselk(-rn + 2, x) - 2*(-rn + 1)*besselk(-rn + 1, x)/x)

    assert check(expand_func(bessely(rn, x)), \
        -bessely(rn - 2, x) + 2*(rn - 1)*bessely(rn - 1, x)/x)
    assert check(expand_func(bessely(-rn, x)), \
        -bessely(-rn + 2, x) + 2*(-rn + 1)*bessely(-rn + 1, x)/x)

    n = Symbol('n', integer=True, positive=True)

    assert expand_func(besseli(n + 2, z)) == \
        besseli(n, z) + (-2*n - 2)*(-2*n*besseli(n, z)/z + besseli(n - 1, z))/z
    assert expand_func(besselj(n + 2, z)) == \
        -besselj(n, z) + (2*n + 2)*(2*n*besselj(n, z)/z - besselj(n - 1, z))/z
    assert expand_func(besselk(n + 2, z)) == \
        besselk(n, z) + (2*n + 2)*(2*n*besselk(n, z)/z + besselk(n - 1, z))/z
    assert expand_func(bessely(n + 2, z)) == \
        -bessely(n, z) + (2*n + 2)*(2*n*bessely(n, z)/z - bessely(n - 1, z))/z

    assert expand_func(besseli(n + S(1)/2, z).rewrite(jn)) == \
        sqrt(2)*sqrt(z)*exp(-I*pi*(n + S(1)/2)/2)* \
        exp_polar(I*pi/4)*jn(n, z*exp_polar(I*pi/2))/sqrt(pi)
    assert expand_func(besselj(n + S(1)/2, z).rewrite(jn)) == \
        sqrt(2)*sqrt(z)*jn(n, z)/sqrt(pi)


def test_fn():
    x, z = symbols("x z")
    assert fn(1, z) == 1/z**2
    assert fn(2, z) == -1/z + 3/z**3
    assert fn(3, z) == -6/z**2 + 15/z**4
    assert fn(4, z) == 1/z - 45/z**3 + 105/z**5


def mjn(n, z):
    return expand_func(jn(n, z))


def myn(n, z):
    return expand_func(yn(n, z))


def test_jn():
    z = symbols("z")
    assert mjn(0, z) == sin(z)/z
    assert mjn(1, z) == sin(z)/z**2 - cos(z)/z
    assert mjn(2, z) == (3/z**3 - 1/z)*sin(z) - (3/z**2) * cos(z)
    assert mjn(3, z) == (15/z**4 - 6/z**2)*sin(z) + (1/z - 15/z**3)*cos(z)
    assert mjn(4, z) == (1/z + 105/z**5 - 45/z**3)*sin(z) + \
        (-105/z**4 + 10/z**2)*cos(z)
    assert mjn(5, z) == (945/z**6 - 420/z**4 + 15/z**2)*sin(z) + \
        (-1/z - 945/z**5 + 105/z**3)*cos(z)
    assert mjn(6, z) == (-1/z + 10395/z**7 - 4725/z**5 + 210/z**3)*sin(z) + \
        (-10395/z**6 + 1260/z**4 - 21/z**2)*cos(z)

    assert expand_func(jn(n, z)) == jn(n, z)


def test_yn():
    z = symbols("z")
    assert myn(0, z) == -cos(z)/z
    assert myn(1, z) == -cos(z)/z**2 - sin(z)/z
    assert myn(2, z) == -((3/z**3 - 1/z)*cos(z) + (3/z**2)*sin(z))
    assert expand_func(yn(n, z)) == yn(n, z)


def test_sympify_yn():
    assert S(15) in myn(3, pi).atoms()
    assert myn(3, pi) == 15/pi**4 - 6/pi**2


def eq(a, b, tol=1e-6):
    for x, y in zip(a, b):
        if not (abs(x - y) < tol):
            return False
    return True


def test_jn_zeros():
    assert eq(jn_zeros(0, 4), [3.141592, 6.283185, 9.424777, 12.566370])
    assert eq(jn_zeros(1, 4), [4.493409, 7.725251, 10.904121, 14.066193])
    assert eq(jn_zeros(2, 4), [5.763459, 9.095011, 12.322940, 15.514603])
    assert eq(jn_zeros(3, 4), [6.987932, 10.417118, 13.698023, 16.923621])
    assert eq(jn_zeros(4, 4), [8.182561, 11.704907, 15.039664, 18.301255])


def test_bessel_eval():
    from sympy import I, Symbol
    n, m, k = Symbol('n', integer=True), Symbol('m'), Symbol('k', integer=True, zero=False)

    for f in [besselj, besseli]:
        assert f(0, 0) == S.One
        assert f(2.1, 0) == S.Zero
        assert f(-3, 0) == S.Zero
        assert f(-10.2, 0) == S.ComplexInfinity
        assert f(1 + 3*I, 0) == S.Zero
        assert f(-3 + I, 0) == S.ComplexInfinity
        assert f(-2*I, 0) == S.NaN
        assert f(n, 0) != S.One and f(n, 0) != S.Zero
        assert f(m, 0) != S.One and f(m, 0) != S.Zero
        assert f(k, 0) == S.Zero

    assert bessely(0, 0) == S.NegativeInfinity
    assert besselk(0, 0) == S.Infinity
    for f in [bessely, besselk]:
        assert f(1 + I, 0) == S.ComplexInfinity
        assert f(I, 0) == S.NaN

    for f in [besselj, bessely]:
        assert f(m, S.Infinity) == S.Zero
        assert f(m, S.NegativeInfinity) == S.Zero

    for f in [besseli, besselk]:
        assert f(m, I*S.Infinity) == S.Zero
        assert f(m, I*S.NegativeInfinity) == S.Zero

    for f in [besseli, besselk]:
        assert f(-4, z) == f(4, z)
        assert f(-3, z) == f(3, z)
        assert f(-n, z) == f(n, z)
        assert f(-m, z) != f(m, z)

    for f in [besselj, bessely]:
        assert f(-4, z) == f(4, z)
        assert f(-3, z) == -f(3, z)
        assert f(-n, z) == (-1)**n*f(n, z)
        assert f(-m, z) != (-1)**m*f(m, z)

    for f in [besselj, besseli]:
        assert f(m, -z) == (-z)**m*z**(-m)*f(m, z)

    assert besseli(2, -z) == besseli(2, z)
    assert besseli(3, -z) == -besseli(3, z)

    assert besselj(0, -z) == besselj(0, z)
    assert besselj(1, -z) == -besselj(1, z)

    assert besseli(0, I*z) == besselj(0, z)
    assert besseli(1, I*z) == I*besselj(1, z)
    assert besselj(3, I*z) == -I*besseli(3, z)


def test_conjugate():
    from sympy import conjugate, I, Symbol
    n, z, x = Symbol('n'), Symbol('z', real=False), Symbol('x', real=True)
    y, t = Symbol('y', real=True, positive=True), Symbol('t', negative=True)

    for f in [besseli, besselj, besselk, bessely, jn, yn, hankel1, hankel2]:
        assert f(n, -1).conjugate() != f(conjugate(n), -1)
        assert f(n, x).conjugate() != f(conjugate(n), x)
        assert f(n, t).conjugate() != f(conjugate(n), t)

    rz = randcplx(b=0.5)

    for f in [besseli, besselj, besselk, bessely, jn, yn]:
        assert f(n, 1 + I).conjugate() == f(conjugate(n), 1 - I)
        assert f(n, 0).conjugate() == f(conjugate(n), 0)
        assert f(n, 1).conjugate() == f(conjugate(n), 1)
        assert f(n, z).conjugate() == f(conjugate(n), conjugate(z))
        assert f(n, y).conjugate() == f(conjugate(n), y)
        assert tn(f(n, rz).conjugate(), f(conjugate(n), conjugate(rz)))

    assert hankel1(n, 1 + I).conjugate() == hankel2(conjugate(n), 1 - I)
    assert hankel1(n, 0).conjugate() == hankel2(conjugate(n), 0)
    assert hankel1(n, 1).conjugate() == hankel2(conjugate(n), 1)
    assert hankel1(n, y).conjugate() == hankel2(conjugate(n), y)
    assert hankel1(n, z).conjugate() == hankel2(conjugate(n), conjugate(z))
    assert tn(hankel1(n, rz).conjugate(), hankel2(conjugate(n), conjugate(rz)))

    assert hankel2(n, 1 + I).conjugate() == hankel1(conjugate(n), 1 - I)
    assert hankel2(n, 0).conjugate() == hankel1(conjugate(n), 0)
    assert hankel2(n, 1).conjugate() == hankel1(conjugate(n), 1)
    assert hankel2(n, y).conjugate() == hankel1(conjugate(n), y)
    assert hankel2(n, z).conjugate() == hankel1(conjugate(n), conjugate(z))
    assert tn(hankel2(n, rz).conjugate(), hankel1(conjugate(n), conjugate(rz)))


def test_branching():
    from sympy import exp_polar, polar_lift, Symbol, I, exp
    assert besselj(polar_lift(k), x) == besselj(k, x)
    assert besseli(polar_lift(k), x) == besseli(k, x)

    n = Symbol('n', integer=True)
    assert besselj(n, exp_polar(2*pi*I)*x) == besselj(n, x)
    assert besselj(n, polar_lift(x)) == besselj(n, x)
    assert besseli(n, exp_polar(2*pi*I)*x) == besseli(n, x)
    assert besseli(n, polar_lift(x)) == besseli(n, x)

    def tn(func, s):
        from random import uniform
        c = uniform(1, 5)
        expr = func(s, c*exp_polar(I*pi)) - func(s, c*exp_polar(-I*pi))
        eps = 1e-15
        expr2 = func(s + eps, -c + eps*I) - func(s + eps, -c - eps*I)
        return abs(expr.n() - expr2.n()).n() < 1e-10

    nu = Symbol('nu')
    assert besselj(nu, exp_polar(2*pi*I)*x) == exp(2*pi*I*nu)*besselj(nu, x)
    assert besseli(nu, exp_polar(2*pi*I)*x) == exp(2*pi*I*nu)*besseli(nu, x)
    assert tn(besselj, 2)
    assert tn(besselj, pi)
    assert tn(besselj, I)
    assert tn(besseli, 2)
    assert tn(besseli, pi)
    assert tn(besseli, I)
