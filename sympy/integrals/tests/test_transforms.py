from sympy.integrals.transforms import (mellin_transform,
    inverse_mellin_transform, laplace_transform,
    inverse_laplace_transform, fourier_transform, inverse_fourier_transform,
    sine_transform, inverse_sine_transform,
    cosine_transform, inverse_cosine_transform,
    hankel_transform, inverse_hankel_transform,
    LaplaceTransform, FourierTransform, SineTransform, CosineTransform,
    InverseLaplaceTransform, InverseFourierTransform,
    InverseSineTransform, InverseCosineTransform, IntegralTransformError)
from sympy.core.function import (Function, expand_mul)
from sympy.core import EulerGamma, Subs, Derivative, diff
from sympy.core.numbers import (I, Rational, oo, pi)
from sympy.core.relational import Eq
from sympy.core.singleton import S
from sympy.core.symbol import (Symbol, symbols)
from sympy.functions.combinatorial.factorials import factorial
from sympy.functions.elementary.complexes import (Abs, re, unpolarify)
from sympy.functions.elementary.exponential import (exp, exp_polar, log)
from sympy.functions.elementary.hyperbolic import (cosh, sinh, coth, asinh)
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import (atan, atan2, cos, sin, tan)
from sympy.functions.special.bessel import (besseli, besselj, besselk, bessely)
from sympy.functions.special.delta_functions import Heaviside
from sympy.functions.special.error_functions import (erf, erfc, expint, Ei)
from sympy.functions.special.gamma_functions import gamma
from sympy.functions.special.hyper import meijerg
from sympy.simplify.gammasimp import gammasimp
from sympy.simplify.hyperexpand import hyperexpand
from sympy.simplify.trigsimp import trigsimp
from sympy.testing.pytest import XFAIL, slow, skip, raises, warns_deprecated_sympy
from sympy.matrices import Matrix, eye
from sympy.abc import x, s, a, b, c, d
nu, beta, rho = symbols('nu beta rho')


def test_undefined_function():
    from sympy.integrals.transforms import MellinTransform
    f = Function('f')
    assert mellin_transform(f(x), x, s) == MellinTransform(f(x), x, s)
    assert mellin_transform(f(x) + exp(-x), x, s) == \
        (MellinTransform(f(x), x, s) + gamma(s + 1)/s, (0, oo), True)


def test_free_symbols():
    f = Function('f')
    assert mellin_transform(f(x), x, s).free_symbols == {s}
    assert mellin_transform(f(x)*a, x, s).free_symbols == {s, a}


def test_as_integral():
    from sympy.integrals.integrals import Integral
    f = Function('f')
    assert mellin_transform(f(x), x, s).rewrite('Integral') == \
        Integral(x**(s - 1)*f(x), (x, 0, oo))
    assert fourier_transform(f(x), x, s).rewrite('Integral') == \
        Integral(f(x)*exp(-2*I*pi*s*x), (x, -oo, oo))
    assert laplace_transform(f(x), x, s).rewrite('Integral') == \
        Integral(f(x)*exp(-s*x), (x, 0, oo))
    assert str(2*pi*I*inverse_mellin_transform(f(s), s, x, (a, b)).rewrite('Integral')) \
        == "Integral(f(s)/x**s, (s, _c - oo*I, _c + oo*I))"
    assert str(2*pi*I*inverse_laplace_transform(f(s), s, x).rewrite('Integral')) == \
        "Integral(f(s)*exp(s*x), (s, _c - oo*I, _c + oo*I))"
    assert inverse_fourier_transform(f(s), s, x).rewrite('Integral') == \
        Integral(f(s)*exp(2*I*pi*s*x), (s, -oo, oo))

# NOTE this is stuck in risch because meijerint cannot handle it


@slow
@XFAIL
def test_mellin_transform_fail():
    skip("Risch takes forever.")

    MT = mellin_transform

    bpos = symbols('b', positive=True)
    # bneg = symbols('b', negative=True)

    expr = (sqrt(x + b**2) + b)**a/sqrt(x + b**2)
    # TODO does not work with bneg, argument wrong. Needs changes to matching.
    assert MT(expr.subs(b, -bpos), x, s) == \
        ((-1)**(a + 1)*2**(a + 2*s)*bpos**(a + 2*s - 1)*gamma(a + s)
         *gamma(1 - a - 2*s)/gamma(1 - s),
            (-re(a), -re(a)/2 + S.Half), True)

    expr = (sqrt(x + b**2) + b)**a
    assert MT(expr.subs(b, -bpos), x, s) == \
        (
            2**(a + 2*s)*a*bpos**(a + 2*s)*gamma(-a - 2*
                   s)*gamma(a + s)/gamma(-s + 1),
            (-re(a), -re(a)/2), True)

    # Test exponent 1:
    assert MT(expr.subs({b: -bpos, a: 1}), x, s) == \
        (-bpos**(2*s + 1)*gamma(s)*gamma(-s - S.Half)/(2*sqrt(pi)),
            (-1, Rational(-1, 2)), True)


def test_mellin_transform():
    from sympy.functions.elementary.miscellaneous import (Max, Min)
    MT = mellin_transform

    bpos = symbols('b', positive=True)

    # 8.4.2
    assert MT(x**nu*Heaviside(x - 1), x, s) == \
        (-1/(nu + s), (-oo, -re(nu)), True)
    assert MT(x**nu*Heaviside(1 - x), x, s) == \
        (1/(nu + s), (-re(nu), oo), True)

    assert MT((1 - x)**(beta - 1)*Heaviside(1 - x), x, s) == \
        (gamma(beta)*gamma(s)/gamma(beta + s), (0, oo), re(beta) > 0)
    assert MT((x - 1)**(beta - 1)*Heaviside(x - 1), x, s) == \
        (gamma(beta)*gamma(1 - beta - s)/gamma(1 - s),
            (-oo, 1 - re(beta)), re(beta) > 0)

    assert MT((1 + x)**(-rho), x, s) == \
        (gamma(s)*gamma(rho - s)/gamma(rho), (0, re(rho)), True)

    assert MT(abs(1 - x)**(-rho), x, s) == (
        2*sin(pi*rho/2)*gamma(1 - rho)*
        cos(pi*(s - rho/2))*gamma(s)*gamma(rho-s)/pi,
        (0, re(rho)), re(rho) < 1)
    mt = MT((1 - x)**(beta - 1)*Heaviside(1 - x)
            + a*(x - 1)**(beta - 1)*Heaviside(x - 1), x, s)
    assert mt[1], mt[2] == ((0, -re(beta) + 1), re(beta) > 0)

    assert MT((x**a - b**a)/(x - b), x, s)[0] == \
        pi*b**(a + s - 1)*sin(pi*a)/(sin(pi*s)*sin(pi*(a + s)))
    assert MT((x**a - bpos**a)/(x - bpos), x, s) == \
        (pi*bpos**(a + s - 1)*sin(pi*a)/(sin(pi*s)*sin(pi*(a + s))),
            (Max(0, -re(a)), Min(1, 1 - re(a))), True)

    expr = (sqrt(x + b**2) + b)**a
    assert MT(expr.subs(b, bpos), x, s) == \
        (-a*(2*bpos)**(a + 2*s)*gamma(s)*gamma(-a - 2*s)/gamma(-a - s + 1),
            (0, -re(a)/2), True)

    expr = (sqrt(x + b**2) + b)**a/sqrt(x + b**2)
    assert MT(expr.subs(b, bpos), x, s) == \
        (2**(a + 2*s)*bpos**(a + 2*s - 1)*gamma(s)
                                         *gamma(1 - a - 2*s)/gamma(1 - a - s),
            (0, -re(a)/2 + S.Half), True)

    # 8.4.2
    assert MT(exp(-x), x, s) == (gamma(s), (0, oo), True)
    assert MT(exp(-1/x), x, s) == (gamma(-s), (-oo, 0), True)

    # 8.4.5
    assert MT(log(x)**4*Heaviside(1 - x), x, s) == (24/s**5, (0, oo), True)
    assert MT(log(x)**3*Heaviside(x - 1), x, s) == (6/s**4, (-oo, 0), True)
    assert MT(log(x + 1), x, s) == (pi/(s*sin(pi*s)), (-1, 0), True)
    assert MT(log(1/x + 1), x, s) == (pi/(s*sin(pi*s)), (0, 1), True)
    assert MT(log(abs(1 - x)), x, s) == (pi/(s*tan(pi*s)), (-1, 0), True)
    assert MT(log(abs(1 - 1/x)), x, s) == (pi/(s*tan(pi*s)), (0, 1), True)

    # 8.4.14
    assert MT(erf(sqrt(x)), x, s) == \
        (-gamma(s + S.Half)/(sqrt(pi)*s), (Rational(-1, 2), 0), True)


def test_mellin_transform2():
    MT = mellin_transform
    # TODO we cannot currently do these (needs summation of 3F2(-1))
    #      this also implies that they cannot be written as a single g-function
    #      (although this is possible)
    mt = MT(log(x)/(x + 1), x, s)
    assert mt[1:] == ((0, 1), True)
    assert not hyperexpand(mt[0], allow_hyper=True).has(meijerg)
    mt = MT(log(x)**2/(x + 1), x, s)
    assert mt[1:] == ((0, 1), True)
    assert not hyperexpand(mt[0], allow_hyper=True).has(meijerg)
    mt = MT(log(x)/(x + 1)**2, x, s)
    assert mt[1:] == ((0, 2), True)
    assert not hyperexpand(mt[0], allow_hyper=True).has(meijerg)


@slow
def test_mellin_transform_bessel():
    from sympy.functions.elementary.miscellaneous import Max
    MT = mellin_transform

    # 8.4.19
    assert MT(besselj(a, 2*sqrt(x)), x, s) == \
        (gamma(a/2 + s)/gamma(a/2 - s + 1), (-re(a)/2, Rational(3, 4)), True)
    assert MT(sin(sqrt(x))*besselj(a, sqrt(x)), x, s) == \
        (2**a*gamma(-2*s + S.Half)*gamma(a/2 + s + S.Half)/(
        gamma(-a/2 - s + 1)*gamma(a - 2*s + 1)), (
        -re(a)/2 - S.Half, Rational(1, 4)), True)
    assert MT(cos(sqrt(x))*besselj(a, sqrt(x)), x, s) == \
        (2**a*gamma(a/2 + s)*gamma(-2*s + S.Half)/(
        gamma(-a/2 - s + S.Half)*gamma(a - 2*s + 1)), (
        -re(a)/2, Rational(1, 4)), True)
    assert MT(besselj(a, sqrt(x))**2, x, s) == \
        (gamma(a + s)*gamma(S.Half - s)
         / (sqrt(pi)*gamma(1 - s)*gamma(1 + a - s)),
            (-re(a), S.Half), True)
    assert MT(besselj(a, sqrt(x))*besselj(-a, sqrt(x)), x, s) == \
        (gamma(s)*gamma(S.Half - s)
         / (sqrt(pi)*gamma(1 - a - s)*gamma(1 + a - s)),
            (0, S.Half), True)
    # NOTE: prudnikov gives the strip below as (1/2 - re(a), 1). As far as
    #       I can see this is wrong (since besselj(z) ~ 1/sqrt(z) for z large)
    assert MT(besselj(a - 1, sqrt(x))*besselj(a, sqrt(x)), x, s) == \
        (gamma(1 - s)*gamma(a + s - S.Half)
         / (sqrt(pi)*gamma(Rational(3, 2) - s)*gamma(a - s + S.Half)),
            (S.Half - re(a), S.Half), True)
    assert MT(besselj(a, sqrt(x))*besselj(b, sqrt(x)), x, s) == \
        (4**s*gamma(1 - 2*s)*gamma((a + b)/2 + s)
         / (gamma(1 - s + (b - a)/2)*gamma(1 - s + (a - b)/2)
            *gamma( 1 - s + (a + b)/2)),
            (-(re(a) + re(b))/2, S.Half), True)
    assert MT(besselj(a, sqrt(x))**2 + besselj(-a, sqrt(x))**2, x, s)[1:] == \
        ((Max(re(a), -re(a)), S.Half), True)

    # Section 8.4.20
    assert MT(bessely(a, 2*sqrt(x)), x, s) == \
        (-cos(pi*(a/2 - s))*gamma(s - a/2)*gamma(s + a/2)/pi,
            (Max(-re(a)/2, re(a)/2), Rational(3, 4)), True)
    assert MT(sin(sqrt(x))*bessely(a, sqrt(x)), x, s) == \
        (-4**s*sin(pi*(a/2 - s))*gamma(S.Half - 2*s)
         * gamma((1 - a)/2 + s)*gamma((1 + a)/2 + s)
         / (sqrt(pi)*gamma(1 - s - a/2)*gamma(1 - s + a/2)),
            (Max(-(re(a) + 1)/2, (re(a) - 1)/2), Rational(1, 4)), True)
    assert MT(cos(sqrt(x))*bessely(a, sqrt(x)), x, s) == \
        (-4**s*cos(pi*(a/2 - s))*gamma(s - a/2)*gamma(s + a/2)*gamma(S.Half - 2*s)
         / (sqrt(pi)*gamma(S.Half - s - a/2)*gamma(S.Half - s + a/2)),
            (Max(-re(a)/2, re(a)/2), Rational(1, 4)), True)
    assert MT(besselj(a, sqrt(x))*bessely(a, sqrt(x)), x, s) == \
        (-cos(pi*s)*gamma(s)*gamma(a + s)*gamma(S.Half - s)
         / (pi**S('3/2')*gamma(1 + a - s)),
            (Max(-re(a), 0), S.Half), True)
    assert MT(besselj(a, sqrt(x))*bessely(b, sqrt(x)), x, s) == \
        (-4**s*cos(pi*(a/2 - b/2 + s))*gamma(1 - 2*s)
         * gamma(a/2 - b/2 + s)*gamma(a/2 + b/2 + s)
         / (pi*gamma(a/2 - b/2 - s + 1)*gamma(a/2 + b/2 - s + 1)),
            (Max((-re(a) + re(b))/2, (-re(a) - re(b))/2), S.Half), True)
    # NOTE bessely(a, sqrt(x))**2 and bessely(a, sqrt(x))*bessely(b, sqrt(x))
    # are a mess (no matter what way you look at it ...)
    assert MT(bessely(a, sqrt(x))**2, x, s)[1:] == \
             ((Max(-re(a), 0, re(a)), S.Half), True)

    # Section 8.4.22
    # TODO we can't do any of these (delicate cancellation)

    # Section 8.4.23
    assert MT(besselk(a, 2*sqrt(x)), x, s) == \
        (gamma(
         s - a/2)*gamma(s + a/2)/2, (Max(-re(a)/2, re(a)/2), oo), True)
    assert MT(besselj(a, 2*sqrt(2*sqrt(x)))*besselk(
        a, 2*sqrt(2*sqrt(x))), x, s) == (4**(-s)*gamma(2*s)*
        gamma(a/2 + s)/(2*gamma(a/2 - s + 1)), (Max(0, -re(a)/2), oo), True)
    # TODO bessely(a, x)*besselk(a, x) is a mess
    assert MT(besseli(a, sqrt(x))*besselk(a, sqrt(x)), x, s) == \
        (gamma(s)*gamma(
        a + s)*gamma(-s + S.Half)/(2*sqrt(pi)*gamma(a - s + 1)),
        (Max(-re(a), 0), S.Half), True)
    assert MT(besseli(b, sqrt(x))*besselk(a, sqrt(x)), x, s) == \
        (2**(2*s - 1)*gamma(-2*s + 1)*gamma(-a/2 + b/2 + s)* \
        gamma(a/2 + b/2 + s)/(gamma(-a/2 + b/2 - s + 1)* \
        gamma(a/2 + b/2 - s + 1)), (Max(-re(a)/2 - re(b)/2, \
        re(a)/2 - re(b)/2), S.Half), True)

    # TODO products of besselk are a mess

    mt = MT(exp(-x/2)*besselk(a, x/2), x, s)
    mt0 = gammasimp(trigsimp(gammasimp(mt[0].expand(func=True))))
    assert mt0 == 2*pi**Rational(3, 2)*cos(pi*s)*gamma(S.Half - s)/(
        (cos(2*pi*a) - cos(2*pi*s))*gamma(-a - s + 1)*gamma(a - s + 1))
    assert mt[1:] == ((Max(-re(a), re(a)), oo), True)
    # TODO exp(x/2)*besselk(a, x/2) [etc] cannot currently be done
    # TODO various strange products of special orders


@slow
def test_expint():
    from sympy.functions.elementary.miscellaneous import Max
    from sympy.functions.special.error_functions import (Ci, E1, Ei, Si)
    from sympy.functions.special.zeta_functions import lerchphi
    from sympy.simplify.simplify import simplify
    aneg = Symbol('a', negative=True)
    u = Symbol('u', polar=True)

    assert mellin_transform(E1(x), x, s) == (gamma(s)/s, (0, oo), True)
    assert inverse_mellin_transform(gamma(s)/s, s, x,
              (0, oo)).rewrite(expint).expand() == E1(x)
    assert mellin_transform(expint(a, x), x, s) == \
        (gamma(s)/(a + s - 1), (Max(1 - re(a), 0), oo), True)
    # XXX IMT has hickups with complicated strips ...
    assert simplify(unpolarify(
                    inverse_mellin_transform(gamma(s)/(aneg + s - 1), s, x,
                  (1 - aneg, oo)).rewrite(expint).expand(func=True))) == \
        expint(aneg, x)

    assert mellin_transform(Si(x), x, s) == \
        (-2**s*sqrt(pi)*gamma(s/2 + S.Half)/(
        2*s*gamma(-s/2 + 1)), (-1, 0), True)
    assert inverse_mellin_transform(-2**s*sqrt(pi)*gamma((s + 1)/2)
                                    /(2*s*gamma(-s/2 + 1)), s, x, (-1, 0)) \
        == Si(x)

    assert mellin_transform(Ci(sqrt(x)), x, s) == \
        (-2**(2*s - 1)*sqrt(pi)*gamma(s)/(s*gamma(-s + S.Half)), (0, 1), True)
    assert inverse_mellin_transform(
        -4**s*sqrt(pi)*gamma(s)/(2*s*gamma(-s + S.Half)),
        s, u, (0, 1)).expand() == Ci(sqrt(u))

    # TODO LT of Si, Shi, Chi is a mess ...
    assert laplace_transform(Ci(x), x, s) == (-log(1 + s**2)/2/s, 0, True)
    assert laplace_transform(expint(a, x), x, s) == \
        (lerchphi(s*exp_polar(I*pi), 1, a), 0, re(a) > S.Zero)
    assert laplace_transform(expint(1, x), x, s) == (log(s + 1)/s, 0, True)
    assert laplace_transform(expint(2, x), x, s) == \
        ((s - log(s + 1))/s**2, 0, True)

    assert inverse_laplace_transform(-log(1 + s**2)/2/s, s, u).expand() == \
        Heaviside(u)*Ci(u)
    assert inverse_laplace_transform(log(s + 1)/s, s, x).rewrite(expint) == \
        Heaviside(x)*E1(x)
    assert inverse_laplace_transform((s - log(s + 1))/s**2, s,
                x).rewrite(expint).expand() == \
        (expint(2, x)*Heaviside(x)).rewrite(Ei).rewrite(expint).expand()


@slow
def test_inverse_mellin_transform():
    from sympy.core.function import expand
    from sympy.functions.elementary.miscellaneous import (Max, Min)
    from sympy.functions.elementary.trigonometric import cot
    from sympy.simplify.powsimp import powsimp
    from sympy.simplify.simplify import simplify
    IMT = inverse_mellin_transform

    assert IMT(gamma(s), s, x, (0, oo)) == exp(-x)
    assert IMT(gamma(-s), s, x, (-oo, 0)) == exp(-1/x)
    assert simplify(IMT(s/(2*s**2 - 2), s, x, (2, oo))) == \
        (x**2 + 1)*Heaviside(1 - x)/(4*x)

    # test passing "None"
    assert IMT(1/(s**2 - 1), s, x, (-1, None)) == \
        -x*Heaviside(-x + 1)/2 - Heaviside(x - 1)/(2*x)
    assert IMT(1/(s**2 - 1), s, x, (None, 1)) == \
        -x*Heaviside(-x + 1)/2 - Heaviside(x - 1)/(2*x)

    # test expansion of sums
    assert IMT(gamma(s) + gamma(s - 1), s, x, (1, oo)) == (x + 1)*exp(-x)/x

    # test factorisation of polys
    r = symbols('r', real=True)
    assert IMT(1/(s**2 + 1), s, exp(-x), (None, oo)
              ).subs(x, r).rewrite(sin).simplify() \
        == sin(r)*Heaviside(1 - exp(-r))

    # test multiplicative substitution
    _a, _b = symbols('a b', positive=True)
    assert IMT(_b**(-s/_a)*factorial(s/_a)/s, s, x, (0, oo)) == exp(-_b*x**_a)
    assert IMT(factorial(_a/_b + s/_b)/(_a + s), s, x, (-_a, oo)) == x**_a*exp(-x**_b)

    def simp_pows(expr):
        return simplify(powsimp(expand_mul(expr, deep=False), force=True)).replace(exp_polar, exp)

    # Now test the inverses of all direct transforms tested above

    # Section 8.4.2
    nu = symbols('nu', real=True)
    assert IMT(-1/(nu + s), s, x, (-oo, None)) == x**nu*Heaviside(x - 1)
    assert IMT(1/(nu + s), s, x, (None, oo)) == x**nu*Heaviside(1 - x)
    assert simp_pows(IMT(gamma(beta)*gamma(s)/gamma(s + beta), s, x, (0, oo))) \
        == (1 - x)**(beta - 1)*Heaviside(1 - x)
    assert simp_pows(IMT(gamma(beta)*gamma(1 - beta - s)/gamma(1 - s),
                         s, x, (-oo, None))) \
        == (x - 1)**(beta - 1)*Heaviside(x - 1)
    assert simp_pows(IMT(gamma(s)*gamma(rho - s)/gamma(rho), s, x, (0, None))) \
        == (1/(x + 1))**rho
    assert simp_pows(IMT(d**c*d**(s - 1)*sin(pi*c)
                         *gamma(s)*gamma(s + c)*gamma(1 - s)*gamma(1 - s - c)/pi,
                         s, x, (Max(-re(c), 0), Min(1 - re(c), 1)))) \
        == (x**c - d**c)/(x - d)

    assert simplify(IMT(1/sqrt(pi)*(-c/2)*gamma(s)*gamma((1 - c)/2 - s)
                        *gamma(-c/2 - s)/gamma(1 - c - s),
                        s, x, (0, -re(c)/2))) == \
        (1 + sqrt(x + 1))**c
    assert simplify(IMT(2**(a + 2*s)*b**(a + 2*s - 1)*gamma(s)*gamma(1 - a - 2*s)
                        /gamma(1 - a - s), s, x, (0, (-re(a) + 1)/2))) == \
        b**(a - 1)*(sqrt(1 + x/b**2) + 1)**(a - 1)*(b**2*sqrt(1 + x/b**2) +
        b**2 + x)/(b**2 + x)
    assert simplify(IMT(-2**(c + 2*s)*c*b**(c + 2*s)*gamma(s)*gamma(-c - 2*s)
                        / gamma(-c - s + 1), s, x, (0, -re(c)/2))) == \
        b**c*(sqrt(1 + x/b**2) + 1)**c

    # Section 8.4.5
    assert IMT(24/s**5, s, x, (0, oo)) == log(x)**4*Heaviside(1 - x)
    assert expand(IMT(6/s**4, s, x, (-oo, 0)), force=True) == \
        log(x)**3*Heaviside(x - 1)
    assert IMT(pi/(s*sin(pi*s)), s, x, (-1, 0)) == log(x + 1)
    assert IMT(pi/(s*sin(pi*s/2)), s, x, (-2, 0)) == log(x**2 + 1)
    assert IMT(pi/(s*sin(2*pi*s)), s, x, (Rational(-1, 2), 0)) == log(sqrt(x) + 1)
    assert IMT(pi/(s*sin(pi*s)), s, x, (0, 1)) == log(1 + 1/x)

    # TODO
    def mysimp(expr):
        from sympy.core.function import expand
        from sympy.simplify.powsimp import powsimp
        from sympy.simplify.simplify import logcombine
        return expand(
            powsimp(logcombine(expr, force=True), force=True, deep=True),
            force=True).replace(exp_polar, exp)

    assert mysimp(mysimp(IMT(pi/(s*tan(pi*s)), s, x, (-1, 0)))) in [
        log(1 - x)*Heaviside(1 - x) + log(x - 1)*Heaviside(x - 1),
        log(x)*Heaviside(x - 1) + log(1 - 1/x)*Heaviside(x - 1) + log(-x +
        1)*Heaviside(-x + 1)]
    # test passing cot
    assert mysimp(IMT(pi*cot(pi*s)/s, s, x, (0, 1))) in [
        log(1/x - 1)*Heaviside(1 - x) + log(1 - 1/x)*Heaviside(x - 1),
        -log(x)*Heaviside(-x + 1) + log(1 - 1/x)*Heaviside(x - 1) + log(-x +
        1)*Heaviside(-x + 1), ]

    # 8.4.14
    assert IMT(-gamma(s + S.Half)/(sqrt(pi)*s), s, x, (Rational(-1, 2), 0)) == \
        erf(sqrt(x))

    # 8.4.19
    assert simplify(IMT(gamma(a/2 + s)/gamma(a/2 - s + 1), s, x, (-re(a)/2, Rational(3, 4)))) \
        == besselj(a, 2*sqrt(x))
    assert simplify(IMT(2**a*gamma(S.Half - 2*s)*gamma(s + (a + 1)/2)
                      / (gamma(1 - s - a/2)*gamma(1 - 2*s + a)),
                      s, x, (-(re(a) + 1)/2, Rational(1, 4)))) == \
        sin(sqrt(x))*besselj(a, sqrt(x))
    assert simplify(IMT(2**a*gamma(a/2 + s)*gamma(S.Half - 2*s)
                      / (gamma(S.Half - s - a/2)*gamma(1 - 2*s + a)),
                      s, x, (-re(a)/2, Rational(1, 4)))) == \
        cos(sqrt(x))*besselj(a, sqrt(x))
    # TODO this comes out as an amazing mess, but simplifies nicely
    assert simplify(IMT(gamma(a + s)*gamma(S.Half - s)
                      / (sqrt(pi)*gamma(1 - s)*gamma(1 + a - s)),
                      s, x, (-re(a), S.Half))) == \
        besselj(a, sqrt(x))**2
    assert simplify(IMT(gamma(s)*gamma(S.Half - s)
                      / (sqrt(pi)*gamma(1 - s - a)*gamma(1 + a - s)),
                      s, x, (0, S.Half))) == \
        besselj(-a, sqrt(x))*besselj(a, sqrt(x))
    assert simplify(IMT(4**s*gamma(-2*s + 1)*gamma(a/2 + b/2 + s)
                      / (gamma(-a/2 + b/2 - s + 1)*gamma(a/2 - b/2 - s + 1)
                         *gamma(a/2 + b/2 - s + 1)),
                      s, x, (-(re(a) + re(b))/2, S.Half))) == \
        besselj(a, sqrt(x))*besselj(b, sqrt(x))

    # Section 8.4.20
    # TODO this can be further simplified!
    assert simplify(IMT(-2**(2*s)*cos(pi*a/2 - pi*b/2 + pi*s)*gamma(-2*s + 1) *
                    gamma(a/2 - b/2 + s)*gamma(a/2 + b/2 + s) /
                    (pi*gamma(a/2 - b/2 - s + 1)*gamma(a/2 + b/2 - s + 1)),
                    s, x,
                    (Max(-re(a)/2 - re(b)/2, -re(a)/2 + re(b)/2), S.Half))) == \
                    besselj(a, sqrt(x))*-(besselj(-b, sqrt(x)) -
                    besselj(b, sqrt(x))*cos(pi*b))/sin(pi*b)
    # TODO more

    # for coverage

    assert IMT(pi/cos(pi*s), s, x, (0, S.Half)) == sqrt(x)/(x + 1)


@slow
def test_laplace_transform():
    from sympy import lowergamma
    from sympy.functions.special.delta_functions import DiracDelta
    from sympy.functions.special.error_functions import (fresnelc, fresnels)
    LT = laplace_transform
    a, b, c, = symbols('a, b, c', positive=True)
    t, w, x = symbols('t, w, x')
    f = Function("f")
    g = Function("g")

    # Test rule-base evaluation according to
    # http://eqworld.ipmnet.ru/en/auxiliary/inttrans/
    # Power-law functions (laplace2.pdf)
    assert LT(a*t+t**2+t**(S(5)/2), t, s) ==\
        (a/s**2 + 2/s**3 + 15*sqrt(pi)/(8*s**(S(7)/2)), 0, True)
    assert LT(b/(t+a), t, s) == (-b*exp(-a*s)*Ei(-a*s), 0, True)
    assert LT(1/sqrt(t+a), t, s) ==\
        (sqrt(pi)*sqrt(1/s)*exp(a*s)*erfc(sqrt(a)*sqrt(s)), 0, True)
    assert LT(sqrt(t)/(t+a), t, s) ==\
        (-pi*sqrt(a)*exp(a*s)*erfc(sqrt(a)*sqrt(s)) + sqrt(pi)*sqrt(1/s),
         0, True)
    assert LT((t+a)**(-S(3)/2), t, s) ==\
        (-2*sqrt(pi)*sqrt(s)*exp(a*s)*erfc(sqrt(a)*sqrt(s)) + 2/sqrt(a),
         0, True)
    assert LT(t**(S(1)/2)*(t+a)**(-1), t, s) ==\
        (-pi*sqrt(a)*exp(a*s)*erfc(sqrt(a)*sqrt(s)) + sqrt(pi)*sqrt(1/s),
         0, True)
    assert LT(1/(a*sqrt(t) + t**(3/2)), t, s) ==\
        (pi*sqrt(a)*exp(a*s)*erfc(sqrt(a)*sqrt(s)), 0, True)
    assert LT((t+a)**b, t, s) ==\
        (s**(-b - 1)*exp(-a*s)*lowergamma(b + 1, a*s), 0, True)
    assert LT(t**5/(t+a), t, s) == (120*a**5*lowergamma(-5, a*s), 0, True)
    # Exponential functions (laplace3.pdf)
    assert LT(exp(t), t, s) == (1/(s - 1), 1, True)
    assert LT(exp(2*t), t, s) == (1/(s - 2), 2, True)
    assert LT(exp(a*t), t, s) == (1/(s - a), a, True)
    assert LT(exp(a*(t-b)), t, s) == (exp(-a*b)/(-a + s), a, True)
    assert LT(t*exp(-a*(t)), t, s) == ((a + s)**(-2), -a, True)
    assert LT(t*exp(-a*(t-b)), t, s) == (exp(a*b)/(a + s)**2, -a, True)
    assert LT(b*t*exp(-a*t), t, s) == (b/(a + s)**2, -a, True)
    assert LT(t**(S(7)/4)*exp(-8*t)/gamma(S(11)/4), t, s) ==\
        ((s + 8)**(-S(11)/4), -8, True)
    assert LT(t**(S(3)/2)*exp(-8*t), t, s) ==\
        (3*sqrt(pi)/(4*(s + 8)**(S(5)/2)), -8, True)
    assert LT(t**a*exp(-a*t), t, s) ==  ((a+s)**(-a-1)*gamma(a+1), -a, True)
    assert LT(b*exp(-a*t**2), t, s) ==\
        (sqrt(pi)*b*exp(s**2/(4*a))*erfc(s/(2*sqrt(a)))/(2*sqrt(a)), 0, True)
    assert LT(exp(-2*t**2), t, s) ==\
        (sqrt(2)*sqrt(pi)*exp(s**2/8)*erfc(sqrt(2)*s/4)/4, 0, True)
    assert LT(b*exp(2*t**2), t, s) == b*LaplaceTransform(exp(2*t**2), t, s)
    assert LT(t*exp(-a*t**2), t, s) ==\
        (1/(2*a) - s*erfc(s/(2*sqrt(a)))/(4*sqrt(pi)*a**(S(3)/2)), 0, True)
    assert LT(exp(-a/t), t, s) ==\
        (2*sqrt(a)*sqrt(1/s)*besselk(1, 2*sqrt(a)*sqrt(s)), 0, True)
    assert LT(sqrt(t)*exp(-a/t), t, s) ==\
        (sqrt(pi)*(2*sqrt(a)*sqrt(s) + 1)*sqrt(s**(-3))*exp(-2*sqrt(a)*\
                                                    sqrt(s))/2, 0, True)
    assert LT(exp(-a/t)/sqrt(t), t, s) ==\
        (sqrt(pi)*sqrt(1/s)*exp(-2*sqrt(a)*sqrt(s)), 0, True)
    assert LT( exp(-a/t)/(t*sqrt(t)), t, s) ==\
        (sqrt(pi)*sqrt(1/a)*exp(-2*sqrt(a)*sqrt(s)), 0, True)
    assert LT(exp(-2*sqrt(a*t)), t, s) ==\
        ( 1/s -sqrt(pi)*sqrt(a) * exp(a/s)*erfc(sqrt(a)*sqrt(1/s))/\
         s**(S(3)/2), 0, True)
    assert LT(exp(-2*sqrt(a*t))/sqrt(t), t, s) == (exp(a/s)*erfc(sqrt(a)*\
        sqrt(1/s))*(sqrt(pi)*sqrt(1/s)), 0, True)
    assert LT(t**4*exp(-2/t), t, s) ==\
        (8*sqrt(2)*(1/s)**(S(5)/2)*besselk(5, 2*sqrt(2)*sqrt(s)), 0, True)
    # Hyperbolic functions (laplace4.pdf)
    assert LT(sinh(a*t), t, s) == (a/(-a**2 + s**2), a, True)
    assert LT(b*sinh(a*t)**2, t, s) == (2*a**2*b/(-4*a**2*s**2 + s**3),
                                        2*a, True)
    # The following line confirms that issue #21202 is solved
    assert LT(cosh(2*t), t, s) == (s/(-4 + s**2), 2, True)
    assert LT(cosh(a*t), t, s) == (s/(-a**2 + s**2), a, True)
    assert LT(cosh(a*t)**2, t, s) == ((-2*a**2 + s**2)/(-4*a**2*s**2 + s**3),
                                      2*a, True)
    assert LT(sinh(x + 3), x, s) == (
        (-s + (s + 1)*exp(6) + 1)*exp(-3)/(s - 1)/(s + 1)/2, 0, Abs(s) > 1)
    # The following line replaces the old test test_issue_7173()
    assert LT(sinh(a*t)*cosh(a*t), t, s) == (a/(-4*a**2 + s**2), 2*a, True)
    assert LT(sinh(a*t)/t, t, s) == (log((a + s)/(-a + s))/2, a, True)
    assert LT(t**(-S(3)/2)*sinh(a*t), t, s) ==\
        (-sqrt(pi)*(sqrt(-a + s) - sqrt(a + s)), a, True)
    assert LT(sinh(2*sqrt(a*t)), t, s) ==\
        (sqrt(pi)*sqrt(a)*exp(a/s)/s**(S(3)/2), 0, True)
    assert LT(sqrt(t)*sinh(2*sqrt(a*t)), t, s) ==\
        (-sqrt(a)/s**2 + sqrt(pi)*(a + s/2)*exp(a/s)*erf(sqrt(a)*\
                                            sqrt(1/s))/s**(S(5)/2), 0, True)
    assert LT(sinh(2*sqrt(a*t))/sqrt(t), t, s) ==\
        (sqrt(pi)*exp(a/s)*erf(sqrt(a)*sqrt(1/s))/sqrt(s), 0, True)
    assert LT(sinh(sqrt(a*t))**2/sqrt(t), t, s) ==\
        (sqrt(pi)*(exp(a/s) - 1)/(2*sqrt(s)), 0, True)
    assert LT(t**(S(3)/7)*cosh(a*t), t, s) ==\
        (((a + s)**(-S(10)/7) + (-a+s)**(-S(10)/7))*gamma(S(10)/7)/2, a, True)
    assert LT(cosh(2*sqrt(a*t)), t, s) ==\
        (sqrt(pi)*sqrt(a)*exp(a/s)*erf(sqrt(a)*sqrt(1/s))/s**(S(3)/2) + 1/s,
         0, True)
    assert LT(sqrt(t)*cosh(2*sqrt(a*t)), t, s) ==\
        (sqrt(pi)*(a + s/2)*exp(a/s)/s**(S(5)/2), 0, True)
    assert LT(cosh(2*sqrt(a*t))/sqrt(t), t, s) ==\
        (sqrt(pi)*exp(a/s)/sqrt(s), 0, True)
    assert LT(cosh(sqrt(a*t))**2/sqrt(t), t, s) ==\
        (sqrt(pi)*(exp(a/s) + 1)/(2*sqrt(s)), 0, True)
    # logarithmic functions (laplace5.pdf)
    assert LT(log(t), t, s) == (-log(s+S.EulerGamma)/s, 0, True)
    assert LT(log(t/a), t, s) == (-log(a*s + S.EulerGamma)/s, 0, True)
    assert LT(log(1+a*t), t, s) == (-exp(s/a)*Ei(-s/a)/s, 0, True)
    assert LT(log(t+a), t, s) == ((log(a) - exp(s/a)*Ei(-s/a)/s)/s, 0, True)
    assert LT(log(t)/sqrt(t), t, s) ==\
        (sqrt(pi)*(-log(s) - 2*log(2) - S.EulerGamma)/sqrt(s), 0, True)
    assert LT(t**(S(5)/2)*log(t), t, s) ==\
        (15*sqrt(pi)*(-log(s)-2*log(2)-S.EulerGamma+S(46)/15)/(8*s**(S(7)/2)),
         0, True)
    assert (LT(t**3*log(t), t, s, noconds=True)-6*(-log(s) - S.EulerGamma\
                                    + S(11)/6)/s**4).simplify() == S.Zero
    assert LT(log(t)**2, t, s) ==\
        (((log(s) + EulerGamma)**2 + pi**2/6)/s, 0, True)
    assert LT(exp(-a*t)*log(t), t, s) ==\
        ((-log(a + s) - S.EulerGamma)/(a + s), -a, True)
    # Trigonometric functions (laplace6.pdf)
    assert LT(sin(a*t), t, s) == (a/(a**2 + s**2), 0, True)
    assert LT(Abs(sin(a*t)), t, s) ==\
        (a*coth(pi*s/(2*a))/(a**2 + s**2), 0, True)
    assert LT(sin(a*t)/t, t, s) == (atan(a/s), 0, True)
    assert LT(sin(a*t)**2/t, t, s) == (log(4*a**2/s**2 + 1)/4, 0, True)
    assert LT(sin(a*t)**2/t**2, t, s) ==\
        (a*atan(2*a/s) - s*log(4*a**2/s**2 + 1)/4, 0, True)
    assert LT(sin(2*sqrt(a*t)), t, s) ==\
        (sqrt(pi)*sqrt(a)*exp(-a/s)/s**(S(3)/2), 0, True)
    assert LT(sin(2*sqrt(a*t))/t, t, s) == (pi*erf(sqrt(a)*sqrt(1/s)), 0, True)
    assert LT(cos(a*t), t, s) == (s/(a**2 + s**2), 0, True)
    assert LT(cos(a*t)**2, t, s) ==\
        ((2*a**2 + s**2)/(s*(4*a**2 + s**2)), 0, True)
    assert LT(sqrt(t)*cos(2*sqrt(a*t)), t, s) ==\
        (sqrt(pi)*(-2*a + s)*exp(-a/s)/(2*s**(S(5)/2)), 0, True)
    assert LT(cos(2*sqrt(a*t))/sqrt(t), t, s) ==\
        (sqrt(pi)*sqrt(1/s)*exp(-a/s), 0, True)
    assert LT(sin(a*t)*sin(b*t), t, s) ==\
        (2*a*b*s/((s**2 + (a - b)**2)*(s**2 + (a + b)**2)), 0, True)
    assert LT(cos(a*t)*sin(b*t), t, s) ==\
        (b*(-a**2 + b**2 + s**2)/((s**2 + (a - b)**2)*(s**2 + (a + b)**2)),
         0, True)
    assert LT(cos(a*t)*cos(b*t), t, s) ==\
        (s*(a**2 + b**2 + s**2)/((s**2 + (a - b)**2)*(s**2 + (a + b)**2)),
         0, True)
    assert LT(c*exp(-b*t)*sin(a*t), t, s) == (a*c/(a**2 + (b + s)**2),
                                              -b, True)
    assert LT(c*exp(-b*t)*cos(a*t), t, s) == ((b + s)*c/(a**2 + (b + s)**2),
                                              -b, True)
    assert LT(cos(x + 3), x, s) == ((s*cos(3) - sin(3))/(s**2 + 1), 0, True)
    # Error functions (laplace7.pdf)
    assert LT(erf(a*t), t, s) == (exp(s**2/(4*a**2))*erfc(s/(2*a))/s, 0, True)
    assert LT(erf(sqrt(a*t)), t, s) == (sqrt(a)/(s*sqrt(a + s)), 0, True)
    assert LT(exp(a*t)*erf(sqrt(a*t)), t, s) ==\
        (sqrt(a)/(sqrt(s)*(-a + s)), a, True)
    assert LT(erf(sqrt(a/t)/2), t, s) == ((1-exp(-sqrt(a)*sqrt(s)))/s, 0, True)
    assert LT(erfc(sqrt(a*t)), t, s) ==\
        ((-sqrt(a) + sqrt(a + s))/(s*sqrt(a + s)), 0, True)
    assert LT(exp(a*t)*erfc(sqrt(a*t)), t, s) ==\
        (1/(sqrt(a)*sqrt(s) + s), 0, True)
    assert LT(erfc(sqrt(a/t)/2), t, s) == (exp(-sqrt(a)*sqrt(s))/s, 0, True)
    # Bessel functions (laplace8.pdf)
    assert LT(besselj(0, a*t), t, s) == (1/sqrt(a**2 + s**2), 0, True)
    assert LT(besselj(1, a*t), t, s) ==\
        (a/(sqrt(a**2 + s**2)*(s + sqrt(a**2 + s**2))), 0, True)
    assert LT(besselj(2, a*t), t, s) ==\
        (a**2/(sqrt(a**2 + s**2)*(s + sqrt(a**2 + s**2))**2), 0, True)
    assert LT(t*besselj(0, a*t), t, s) ==\
        (s/(a**2 + s**2)**(S(3)/2), 0, True)
    assert LT(t*besselj(1, a*t), t, s) ==\
        (a/(a**2 + s**2)**(S(3)/2), 0, True)
    assert LT(t**2*besselj(2, a*t), t, s) ==\
        (3*a**2/(a**2 + s**2)**(S(5)/2), 0, True)
    assert LT(besselj(0, 2*sqrt(a*t)), t, s) == (exp(-a/s)/s, 0, True)
    assert LT(t**(S(3)/2)*besselj(3, 2*sqrt(a*t)), t, s) ==\
        (a**(S(3)/2)*exp(-a/s)/s**4, 0, True)
    assert LT(besselj(0, a*sqrt(t**2+b*t)), t, s) ==\
        (exp(b*s - b*sqrt(a**2 + s**2))/sqrt(a**2 + s**2), 0, True)
    assert LT(besseli(0, a*t), t, s) == (1/sqrt(-a**2 + s**2), a, True)
    assert LT(besseli(1, a*t), t, s) ==\
        (a/(sqrt(-a**2 + s**2)*(s + sqrt(-a**2 + s**2))), a, True)
    assert LT(besseli(2, a*t), t, s) ==\
        (a**2/(sqrt(-a**2 + s**2)*(s + sqrt(-a**2 + s**2))**2), a, True)
    assert LT(t*besseli(0, a*t), t, s) == (s/(-a**2 + s**2)**(S(3)/2), a, True)
    assert LT(t*besseli(1, a*t), t, s) == (a/(-a**2 + s**2)**(S(3)/2), a, True)
    assert LT(t**2*besseli(2, a*t), t, s) ==\
        (3*a**2/(-a**2 + s**2)**(S(5)/2), a, True)
    assert LT(t**(S(3)/2)*besseli(3, 2*sqrt(a*t)), t, s) ==\
        (a**(S(3)/2)*exp(a/s)/s**4, 0, True)
    assert LT(bessely(0, a*t), t, s) ==\
        (-2*asinh(s/a)/(pi*sqrt(a**2 + s**2)), 0, True)
    assert LT(besselk(0, a*t), t, s) ==\
        (log(s + sqrt(-a**2 + s**2))/sqrt(-a**2 + s**2), a, True)
    assert LT(sin(a*t)**8, t, s) ==\
        (40320*a**8/(s*(147456*a**8 + 52480*a**6*s**2 + 4368*a**4*s**4 +\
                        120*a**2*s**6 + s**8)), 0, True)

    # Test general rules and unevaluated forms
    # These all also test whether issue #7219 is solved.
    assert LT(Heaviside(t-1)*cos(t-1), t, s) == (s*exp(-s)/(s**2 + 1), 0, True)
    assert LT(a*f(t), t, w) == a*LaplaceTransform(f(t), t, w)
    assert LT(a*Heaviside(t+1)*f(t+1), t, s) ==\
        a*LaplaceTransform(f(t + 1)*Heaviside(t + 1), t, s)
    assert LT(a*Heaviside(t-1)*f(t-1), t, s) ==\
        a*LaplaceTransform(f(t), t, s)*exp(-s)
    assert LT(b*f(t/a), t, s) == a*b*LaplaceTransform(f(t), t, a*s)
    assert LT(exp(-f(x)*t), t, s) == (1/(s + f(x)), -f(x), True)
    assert LT(exp(-a*t)*f(t), t, s) == LaplaceTransform(f(t), t, a + s)
    assert LT(exp(-a*t)*erfc(sqrt(b/t)/2), t, s) ==\
        (exp(-sqrt(b)*sqrt(a + s))/(a + s), -a, True)
    assert LT(sinh(a*t)*f(t), t, s) ==\
        LaplaceTransform(f(t), t, -a+s)/2 - LaplaceTransform(f(t), t, a+s)/2
    assert LT(sinh(a*t)*t, t, s) ==\
        (-1/(2*(a + s)**2) + 1/(2*(-a + s)**2), a, True)
    assert LT(cosh(a*t)*f(t), t, s) ==\
        LaplaceTransform(f(t), t, -a+s)/2 + LaplaceTransform(f(t), t, a+s)/2
    assert LT(cosh(a*t)*t, t, s) ==\
        (1/(2*(a + s)**2) + 1/(2*(-a + s)**2), a, True)
    assert LT(sin(a*t)*f(t), t, s) ==\
        I*(-LaplaceTransform(f(t), t, -I*a + s) +\
           LaplaceTransform(f(t), t, I*a + s))/2
    assert LT(sin(a*t)*t, t, s) ==\
        (2*a*s/(a**4 + 2*a**2*s**2 + s**4), 0, True)
    assert LT(cos(a*t)*f(t), t, s) ==\
        LaplaceTransform(f(t), t, -I*a + s)/2 +\
        LaplaceTransform(f(t), t, I*a + s)/2
    assert LT(cos(a*t)*t, t, s) ==\
        ((-a**2 + s**2)/(a**4 + 2*a**2*s**2 + s**4), 0, True)
    # The following two lines test whether issues #5813 and #7176 are solved.
    assert LT(diff(f(t), (t, 1)), t, s) == s*LaplaceTransform(f(t), t, s)\
        - f(0)
    assert LT(diff(f(t), (t, 3)), t, s) == s**3*LaplaceTransform(f(t), t, s)\
        - s**2*f(0) - s*Subs(Derivative(f(t), t), t, 0)\
            - Subs(Derivative(f(t), (t, 2)), t, 0)
    assert LT(a*f(b*t)+g(c*t), t, s) == a*LaplaceTransform(f(t), t, s/b)/b +\
        LaplaceTransform(g(t), t, s/c)/c
    assert inverse_laplace_transform(
        f(w), w, t, plane=0) == InverseLaplaceTransform(f(w), w, t, 0)
    assert LT(f(t)*g(t), t, s) == LaplaceTransform(f(t)*g(t), t, s)

    # additional basic tests from wikipedia
    assert LT((t - a)**b*exp(-c*(t - a))*Heaviside(t - a), t, s) == \
        ((s + c)**(-b - 1)*exp(-a*s)*gamma(b + 1), -c, True)
    assert LT((exp(2*t) - 1)*exp(-b - t)*Heaviside(t)/2, t, s, noconds=True) \
        == exp(-b)/(s**2 - 1)

    # DiracDelta function: standard cases
    assert LT(DiracDelta(t), t, s) == (1, 0, True)
    assert LT(DiracDelta(a*t), t, s) == (1/a, 0, True)
    assert LT(DiracDelta(t/42), t, s) == (42, 0, True)
    assert LT(DiracDelta(t+42), t, s) == (0, 0, True)
    assert LT(DiracDelta(t)+DiracDelta(t-42), t, s) == \
        (1 + exp(-42*s), 0, True)
    assert LT(DiracDelta(t)-a*exp(-a*t), t, s) == (s/(a + s), 0, True)
    assert LT(exp(-t)*(DiracDelta(t)+DiracDelta(t-42)), t, s) == \
        (exp(-42*s - 42) + 1, -oo, True)

    # Collection of cases that cannot be fully evaluated and/or would catch
    # some common implementation errors
    assert LT(DiracDelta(t**2), t, s) == LaplaceTransform(DiracDelta(t**2), t, s)
    assert LT(DiracDelta(t**2 - 1), t, s) == (exp(-s)/2, -oo, True)
    assert LT(DiracDelta(t*(1 - t)), t, s) == \
        LaplaceTransform(DiracDelta(-t**2 + t), t, s)
    assert LT((DiracDelta(t) + 1)*(DiracDelta(t - 1) + 1), t, s) == \
        (LaplaceTransform(DiracDelta(t)*DiracDelta(t - 1), t, s) + \
         1 + exp(-s) + 1/s, 0, True)
    assert LT(DiracDelta(2*t-2*exp(a)), t, s) == (exp(-s*exp(a))/2, 0, True)
    assert LT(DiracDelta(-2*t+2*exp(a)), t, s) == (exp(-s*exp(a))/2, 0, True)

    # Heaviside tests
    assert LT(Heaviside(t), t, s) == (1/s, 0, True)
    assert LT(Heaviside(t - a), t, s) == (exp(-a*s)/s, 0, True)
    assert LT(Heaviside(t-1), t, s) == (exp(-s)/s, 0, True)
    assert LT(Heaviside(2*t-4), t, s) == (exp(-2*s)/s, 0, True)
    assert LT(Heaviside(-2*t+4), t, s) == ((1 - exp(-2*s))/s, 0, True)
    assert LT(Heaviside(2*t+4), t, s) == (1/s, 0, True)
    assert LT(Heaviside(-2*t+4), t, s) == ((1 - exp(-2*s))/s, 0, True)

    # Fresnel functions
    assert laplace_transform(fresnels(t), t, s) == \
        ((-sin(s**2/(2*pi))*fresnels(s/pi) + sin(s**2/(2*pi))/2 -
            cos(s**2/(2*pi))*fresnelc(s/pi) + cos(s**2/(2*pi))/2)/s, 0, True)
    assert laplace_transform(fresnelc(t), t, s) == (
        ((2*sin(s**2/(2*pi))*fresnelc(s/pi) - 2*cos(s**2/(2*pi))*fresnels(s/pi)
        + sqrt(2)*cos(s**2/(2*pi) + pi/4))/(2*s), 0, True))

    # Matrix tests
    Mt = Matrix([[exp(t), t*exp(-t)], [t*exp(-t), exp(t)]])
    Ms = Matrix([[    1/(s - 1), (s + 1)**(-2)],
                 [(s + 1)**(-2),     1/(s - 1)]])

    # The default behaviour for Laplace tranform of a Matrix returns a Matrix
    # of Tuples and is deprecated:
    with warns_deprecated_sympy():
        Ms_conds = Matrix([[(1/(s - 1), 1, True), ((s + 1)**(-2),
            -1, True)], [((s + 1)**(-2), -1, True), (1/(s - 1), 1, True)]])
    with warns_deprecated_sympy():
        assert LT(Mt, t, s) == Ms_conds
    # The new behavior is to return a tuple of a Matrix and the convergence
    # conditions for the matrix as a whole:
    assert LT(Mt, t, s, legacy_matrix=False) == (Ms, 1, True)
    # With noconds=True the transformed matrix is returned without conditions
    # either way:
    assert LT(Mt, t, s, noconds=True) == Ms
    assert LT(Mt, t, s, legacy_matrix=False, noconds=True) == Ms

@slow
def test_inverse_laplace_transform():
    from sympy.core.exprtools import factor_terms
    from sympy.functions.special.delta_functions import DiracDelta
    from sympy.simplify.simplify import simplify
    ILT = inverse_laplace_transform
    a, b, c, = symbols('a b c', positive=True)
    t = symbols('t')

    def simp_hyp(expr):
        return factor_terms(expand_mul(expr)).rewrite(sin)

    assert ILT(1, s, t) == DiracDelta(t)
    assert ILT(1/s, s, t) == Heaviside(t)
    assert ILT(a/(a + s), s, t) == a*exp(-a*t)*Heaviside(t)
    assert ILT(s/(a + s), s, t) == -a*exp(-a*t)*Heaviside(t) + DiracDelta(t)
    assert ILT((a + s)**(-2), s, t) == t*exp(-a*t)*Heaviside(t)
    assert ILT((a + s)**(-5), s, t) == t**4*exp(-a*t)*Heaviside(t)/24
    assert ILT(a/(a**2 + s**2), s, t) == sin(a*t)*Heaviside(t)
    assert ILT(s/(s**2 + a**2), s, t) == cos(a*t)*Heaviside(t)
    assert ILT(b/(b**2 + (a + s)**2), s, t) == exp(-a*t)*sin(b*t)*Heaviside(t)
    assert ILT(b*s/(b**2 + (a + s)**2), s, t) +\
        (a*sin(b*t) - b*cos(b*t))*exp(-a*t)*Heaviside(t) == 0
    assert ILT(exp(-a*s)/s, s, t) == Heaviside(-a + t)
    assert ILT(exp(-a*s)/(b + s), s, t) == exp(b*(a - t))*Heaviside(-a + t)
    assert ILT((b + s)/(a**2 + (b + s)**2), s, t) == \
        exp(-b*t)*cos(a*t)*Heaviside(t)
    assert ILT(exp(-a*s)/s**b, s, t) == \
        (-a + t)**(b - 1)*Heaviside(-a + t)/gamma(b)
    assert ILT(exp(-a*s)/sqrt(s**2 + 1), s, t) == \
        Heaviside(-a + t)*besselj(0, a - t)
    assert ILT(1/(s*sqrt(s + 1)), s, t) == Heaviside(t)*erf(sqrt(t))
    assert ILT(1/(s**2*(s**2 + 1)), s, t) == (t - sin(t))*Heaviside(t)
    assert ILT(s**2/(s**2 + 1), s, t) == -sin(t)*Heaviside(t) + DiracDelta(t)
    assert ILT(1 - 1/(s**2 + 1), s, t) == -sin(t)*Heaviside(t) + DiracDelta(t)
    assert ILT(1/s**2, s, t) == t*Heaviside(t)
    assert ILT(1/s**5, s, t) == t**4*Heaviside(t)/24
    assert simp_hyp(ILT(a/(s**2 - a**2), s, t)) == sinh(a*t)*Heaviside(t)
    assert simp_hyp(ILT(s/(s**2 - a**2), s, t)) == cosh(a*t)*Heaviside(t)
    # TODO sinh/cosh shifted come out a mess. also delayed trig is a mess
    # TODO should this simplify further?
    assert ILT(exp(-a*s)/s**b, s, t) == \
        (t - a)**(b - 1)*Heaviside(t - a)/gamma(b)
    assert ILT(exp(-a*s)/sqrt(1 + s**2), s, t) == \
        Heaviside(t - a)*besselj(0, a - t)  # note: besselj(0, x) is even
    # XXX ILT turns these branch factor into trig functions ...
    assert simplify(ILT(a**b*(s + sqrt(s**2 - a**2))**(-b)/sqrt(s**2 - a**2),
                    s, t).rewrite(exp)) == \
        Heaviside(t)*besseli(b, a*t)
    assert ILT(a**b*(s + sqrt(s**2 + a**2))**(-b)/sqrt(s**2 + a**2),
               s, t).rewrite(exp) == \
        Heaviside(t)*besselj(b, a*t)

    assert ILT(1/(s*sqrt(s + 1)), s, t) == Heaviside(t)*erf(sqrt(t))
    # TODO can we make erf(t) work?

    assert ILT(1/(s**2*(s**2 + 1)),s,t) == (t - sin(t))*Heaviside(t)

    assert ILT( (s * eye(2) - Matrix([[1, 0], [0, 2]])).inv(), s, t) ==\
        Matrix([[exp(t)*Heaviside(t), 0], [0, exp(2*t)*Heaviside(t)]])

def test_inverse_laplace_transform_delta():
    from sympy.functions.special.delta_functions import DiracDelta
    ILT = inverse_laplace_transform
    t = symbols('t')
    assert ILT(2, s, t) == 2*DiracDelta(t)
    assert ILT(2*exp(3*s) - 5*exp(-7*s), s, t) == \
        2*DiracDelta(t + 3) - 5*DiracDelta(t - 7)
    a = cos(sin(7)/2)
    assert ILT(a*exp(-3*s), s, t) == a*DiracDelta(t - 3)
    assert ILT(exp(2*s), s, t) == DiracDelta(t + 2)
    r = Symbol('r', real=True)
    assert ILT(exp(r*s), s, t) == DiracDelta(t + r)


def test_inverse_laplace_transform_delta_cond():
    from sympy.functions.elementary.complexes import im
    from sympy.functions.special.delta_functions import DiracDelta
    ILT = inverse_laplace_transform
    t = symbols('t')
    r = Symbol('r', real=True)
    assert ILT(exp(r*s), s, t, noconds=False) == (DiracDelta(t + r), True)
    z = Symbol('z')
    assert ILT(exp(z*s), s, t, noconds=False) == \
        (DiracDelta(t + z), Eq(im(z), 0))
    # inversion does not exist: verify it doesn't evaluate to DiracDelta
    for z in (Symbol('z', extended_real=False),
              Symbol('z', imaginary=True, zero=False)):
        f = ILT(exp(z*s), s, t, noconds=False)
        f = f[0] if isinstance(f, tuple) else f
        assert f.func != DiracDelta
    # issue 15043
    assert ILT(1/s + exp(r*s)/s, s, t, noconds=False) == (
        Heaviside(t) + Heaviside(r + t), True)

def test_fourier_transform():
    from sympy.core.function import (expand, expand_complex, expand_trig)
    from sympy.polys.polytools import factor
    from sympy.simplify.simplify import simplify
    FT = fourier_transform
    IFT = inverse_fourier_transform

    def simp(x):
        return simplify(expand_trig(expand_complex(expand(x))))

    def sinc(x):
        return sin(pi*x)/(pi*x)
    k = symbols('k', real=True)
    f = Function("f")

    # TODO for this to work with real a, need to expand abs(a*x) to abs(a)*abs(x)
    a = symbols('a', positive=True)
    b = symbols('b', positive=True)

    posk = symbols('posk', positive=True)

    # Test unevaluated form
    assert fourier_transform(f(x), x, k) == FourierTransform(f(x), x, k)
    assert inverse_fourier_transform(
        f(k), k, x) == InverseFourierTransform(f(k), k, x)

    # basic examples from wikipedia
    assert simp(FT(Heaviside(1 - abs(2*a*x)), x, k)) == sinc(k/a)/a
    # TODO IFT is a *mess*
    assert simp(FT(Heaviside(1 - abs(a*x))*(1 - abs(a*x)), x, k)) == sinc(k/a)**2/a
    # TODO IFT

    assert factor(FT(exp(-a*x)*Heaviside(x), x, k), extension=I) == \
        1/(a + 2*pi*I*k)
    # NOTE: the ift comes out in pieces
    assert IFT(1/(a + 2*pi*I*x), x, posk,
            noconds=False) == (exp(-a*posk), True)
    assert IFT(1/(a + 2*pi*I*x), x, -posk,
            noconds=False) == (0, True)
    assert IFT(1/(a + 2*pi*I*x), x, symbols('k', negative=True),
            noconds=False) == (0, True)
    # TODO IFT without factoring comes out as meijer g

    assert factor(FT(x*exp(-a*x)*Heaviside(x), x, k), extension=I) == \
        1/(a + 2*pi*I*k)**2
    assert FT(exp(-a*x)*sin(b*x)*Heaviside(x), x, k) == \
        b/(b**2 + (a + 2*I*pi*k)**2)

    assert FT(exp(-a*x**2), x, k) == sqrt(pi)*exp(-pi**2*k**2/a)/sqrt(a)
    assert IFT(sqrt(pi/a)*exp(-(pi*k)**2/a), k, x) == exp(-a*x**2)
    assert FT(exp(-a*abs(x)), x, k) == 2*a/(a**2 + 4*pi**2*k**2)
    # TODO IFT (comes out as meijer G)

    # TODO besselj(n, x), n an integer > 0 actually can be done...

    # TODO are there other common transforms (no distributions!)?


def test_sine_transform():
    t = symbols("t")
    w = symbols("w")
    a = symbols("a")
    f = Function("f")

    # Test unevaluated form
    assert sine_transform(f(t), t, w) == SineTransform(f(t), t, w)
    assert inverse_sine_transform(
        f(w), w, t) == InverseSineTransform(f(w), w, t)

    assert sine_transform(1/sqrt(t), t, w) == 1/sqrt(w)
    assert inverse_sine_transform(1/sqrt(w), w, t) == 1/sqrt(t)

    assert sine_transform((1/sqrt(t))**3, t, w) == 2*sqrt(w)

    assert sine_transform(t**(-a), t, w) == 2**(
        -a + S.Half)*w**(a - 1)*gamma(-a/2 + 1)/gamma((a + 1)/2)
    assert inverse_sine_transform(2**(-a + S(
        1)/2)*w**(a - 1)*gamma(-a/2 + 1)/gamma(a/2 + S.Half), w, t) == t**(-a)

    assert sine_transform(
        exp(-a*t), t, w) == sqrt(2)*w/(sqrt(pi)*(a**2 + w**2))
    assert inverse_sine_transform(
        sqrt(2)*w/(sqrt(pi)*(a**2 + w**2)), w, t) == exp(-a*t)

    assert sine_transform(
        log(t)/t, t, w) == sqrt(2)*sqrt(pi)*-(log(w**2) + 2*EulerGamma)/4

    assert sine_transform(
        t*exp(-a*t**2), t, w) == sqrt(2)*w*exp(-w**2/(4*a))/(4*a**Rational(3, 2))
    assert inverse_sine_transform(
        sqrt(2)*w*exp(-w**2/(4*a))/(4*a**Rational(3, 2)), w, t) == t*exp(-a*t**2)


def test_cosine_transform():
    from sympy.functions.special.error_functions import (Ci, Si)

    t = symbols("t")
    w = symbols("w")
    a = symbols("a")
    f = Function("f")

    # Test unevaluated form
    assert cosine_transform(f(t), t, w) == CosineTransform(f(t), t, w)
    assert inverse_cosine_transform(
        f(w), w, t) == InverseCosineTransform(f(w), w, t)

    assert cosine_transform(1/sqrt(t), t, w) == 1/sqrt(w)
    assert inverse_cosine_transform(1/sqrt(w), w, t) == 1/sqrt(t)

    assert cosine_transform(1/(
        a**2 + t**2), t, w) == sqrt(2)*sqrt(pi)*exp(-a*w)/(2*a)

    assert cosine_transform(t**(
        -a), t, w) == 2**(-a + S.Half)*w**(a - 1)*gamma((-a + 1)/2)/gamma(a/2)
    assert inverse_cosine_transform(2**(-a + S(
        1)/2)*w**(a - 1)*gamma(-a/2 + S.Half)/gamma(a/2), w, t) == t**(-a)

    assert cosine_transform(
        exp(-a*t), t, w) == sqrt(2)*a/(sqrt(pi)*(a**2 + w**2))
    assert inverse_cosine_transform(
        sqrt(2)*a/(sqrt(pi)*(a**2 + w**2)), w, t) == exp(-a*t)

    assert cosine_transform(exp(-a*sqrt(t))*cos(a*sqrt(
        t)), t, w) == a*exp(-a**2/(2*w))/(2*w**Rational(3, 2))

    assert cosine_transform(1/(a + t), t, w) == sqrt(2)*(
        (-2*Si(a*w) + pi)*sin(a*w)/2 - cos(a*w)*Ci(a*w))/sqrt(pi)
    assert inverse_cosine_transform(sqrt(2)*meijerg(((S.Half, 0), ()), (
        (S.Half, 0, 0), (S.Half,)), a**2*w**2/4)/(2*pi), w, t) == 1/(a + t)

    assert cosine_transform(1/sqrt(a**2 + t**2), t, w) == sqrt(2)*meijerg(
        ((S.Half,), ()), ((0, 0), (S.Half,)), a**2*w**2/4)/(2*sqrt(pi))
    assert inverse_cosine_transform(sqrt(2)*meijerg(((S.Half,), ()), ((0, 0), (S.Half,)), a**2*w**2/4)/(2*sqrt(pi)), w, t) == 1/(t*sqrt(a**2/t**2 + 1))


def test_hankel_transform():
    r = Symbol("r")
    k = Symbol("k")
    nu = Symbol("nu")
    m = Symbol("m")
    a = symbols("a")

    assert hankel_transform(1/r, r, k, 0) == 1/k
    assert inverse_hankel_transform(1/k, k, r, 0) == 1/r

    assert hankel_transform(
        1/r**m, r, k, 0) == 2**(-m + 1)*k**(m - 2)*gamma(-m/2 + 1)/gamma(m/2)
    assert inverse_hankel_transform(
        2**(-m + 1)*k**(m - 2)*gamma(-m/2 + 1)/gamma(m/2), k, r, 0) == r**(-m)

    assert hankel_transform(1/r**m, r, k, nu) == (
        2*2**(-m)*k**(m - 2)*gamma(-m/2 + nu/2 + 1)/gamma(m/2 + nu/2))
    assert inverse_hankel_transform(2**(-m + 1)*k**(
        m - 2)*gamma(-m/2 + nu/2 + 1)/gamma(m/2 + nu/2), k, r, nu) == r**(-m)

    assert hankel_transform(r**nu*exp(-a*r), r, k, nu) == \
        2**(nu + 1)*a*k**(-nu - 3)*(a**2/k**2 + 1)**(-nu - S(
                                                     3)/2)*gamma(nu + Rational(3, 2))/sqrt(pi)
    assert inverse_hankel_transform(
        2**(nu + 1)*a*k**(-nu - 3)*(a**2/k**2 + 1)**(-nu - Rational(3, 2))*gamma(
        nu + Rational(3, 2))/sqrt(pi), k, r, nu) == r**nu*exp(-a*r)


def test_issue_7181():
    assert mellin_transform(1/(1 - x), x, s) != None


def test_issue_8882():
    # This is the original test.
    # from sympy import diff, Integral, integrate
    # r = Symbol('r')
    # psi = 1/r*sin(r)*exp(-(a0*r))
    # h = -1/2*diff(psi, r, r) - 1/r*psi
    # f = 4*pi*psi*h*r**2
    # assert integrate(f, (r, -oo, 3), meijerg=True).has(Integral) == True

    # To save time, only the critical part is included.
    F = -a**(-s + 1)*(4 + 1/a**2)**(-s/2)*sqrt(1/a**2)*exp(-s*I*pi)* \
        sin(s*atan(sqrt(1/a**2)/2))*gamma(s)
    raises(IntegralTransformError, lambda:
        inverse_mellin_transform(F, s, x, (-1, oo),
        **{'as_meijerg': True, 'needeval': True}))


def test_issue_8514():
    from sympy.simplify.simplify import simplify
    a, b, c, = symbols('a b c', positive=True)
    t = symbols('t', positive=True)
    ft = simplify(inverse_laplace_transform(1/(a*s**2+b*s+c),s, t))
    assert ft == (I*exp(t*cos(atan2(0, -4*a*c + b**2)/2)*sqrt(Abs(4*a*c -
                  b**2))/a)*sin(t*sin(atan2(0, -4*a*c + b**2)/2)*sqrt(Abs(
                  4*a*c - b**2))/(2*a)) + exp(t*cos(atan2(0, -4*a*c + b**2)
                  /2)*sqrt(Abs(4*a*c - b**2))/a)*cos(t*sin(atan2(0, -4*a*c
                  + b**2)/2)*sqrt(Abs(4*a*c - b**2))/(2*a)) + I*sin(t*sin(
                  atan2(0, -4*a*c + b**2)/2)*sqrt(Abs(4*a*c - b**2))/(2*a))
                  - cos(t*sin(atan2(0, -4*a*c + b**2)/2)*sqrt(Abs(4*a*c -
                  b**2))/(2*a)))*exp(-t*(b + cos(atan2(0, -4*a*c + b**2)/2)
                  *sqrt(Abs(4*a*c - b**2)))/(2*a))/sqrt(-4*a*c + b**2)


def test_issue_12591():
    x, y = symbols("x y", real=True)
    assert fourier_transform(exp(x), x, y) == FourierTransform(exp(x), x, y)


def test_issue_14692():
    b = Symbol('b', negative=True)
    assert laplace_transform(1/(I*x - b), x, s) == \
        (-I*exp(I*b*s)*expint(1, b*s*exp_polar(I*pi/2)), 0, True)
