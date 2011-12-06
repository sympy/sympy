from sympy.integrals.transforms import (mellin_transform,
    inverse_mellin_transform, laplace_transform, inverse_laplace_transform,
    fourier_transform, inverse_fourier_transform)
from sympy import (gamma, exp, oo, Heaviside, symbols, re, factorial, pi,
                   cos, S, And, sin, sqrt, I, log, tan, hyperexpand, meijerg,
                   EulerGamma, erf, besselj, bessely, besseli, besselk,
                   exp_polar, polar_lift, unpolarify)
from sympy.abc import x, s, a, b
nu, beta, rho = symbols('nu beta rho')

def test_undefined_function():
    from sympy import Function, MellinTransform
    f = Function('f')
    assert mellin_transform(f(x), x, s) == MellinTransform(f(x), x, s)
    assert mellin_transform(f(x) + exp(-x), x, s) == \
           (MellinTransform(f(x), x, s) + gamma(s), (0, oo), True)

def test_free_symbols():
    from sympy import Function
    f = Function('f')
    assert mellin_transform(f(x), x, s).free_symbols == set([s])
    assert mellin_transform(f(x)*a, x, s).free_symbols == set([s, a])

def test_as_integral():
    from sympy import Function, Integral
    f = Function('f')
    assert mellin_transform(f(x), x, s).rewrite('Integral') == \
           Integral(x**(s - 1)*f(x), (x, 0, oo))

def test_mellin_transform():
    from sympy import Max, Min
    MT = mellin_transform

    bpos = symbols('b', positive=True)

    # 8.4.2
    assert MT(x**nu*Heaviside(x - 1), x, s) \
           == (1/(-nu - s), (-oo, -re(nu)), True)
    assert MT(x**nu*Heaviside(1 - x), x, s) \
           == (1/(nu + s), (-re(nu), oo), True)

    assert MT((1-x)**(beta - 1)*Heaviside(1-x), x, s) \
           == (gamma(beta)*gamma(s)/gamma(beta + s),
               (0, oo), re(-beta) < 0)
    assert MT((x-1)**(beta - 1)*Heaviside(x-1), x, s) \
           == (gamma(beta)*gamma(1 - beta - s)/gamma(1 - s),
               (-oo, -re(beta) + 1), re(-beta) < 0)

    assert MT((1+x)**(-rho), x, s) == (gamma(s)*gamma(rho-s)/gamma(rho),
                                       (0, re(rho)), True)

    # TODO also the conditions should be simplified
    assert MT(abs(1-x)**(-rho), x, s) == \
        (cos(pi*rho/2 - pi*s)*gamma(s)*gamma(rho-s)/(cos(pi*rho/2)*gamma(rho)),\
         (0, re(rho)), And(re(rho) - 1 < 0, re(rho) < 1))

    mt = MT((1-x)**(beta-1)*Heaviside(1-x)
            + a*(x-1)**(beta-1)*Heaviside(x-1), x, s)
    assert mt[1], mt[2] == ((0, -re(beta) + 1), True)

    assert MT((x**a-b**a)/(x-b), x, s)[0] == \
           pi*b**(a+s-1)*sin(pi*a)/(sin(pi*s)*sin(pi*(a + s)))
    assert MT((x**a-bpos**a)/(x-bpos), x, s) == \
           (pi*bpos**(a+s-1)*sin(pi*a)/(sin(pi*s)*sin(pi*(a + s))),
            (Max(-re(a), 0), Min(1 - re(a), 1)), True)

    expr = (sqrt(x+b**2)+b)**a/sqrt(x+b**2)
    assert MT(expr.subs(b, bpos), x, s) == \
           (2**(a + 2*s)*bpos**(a + 2*s - 1)*gamma(s) \
                                         *gamma(1 - a - 2*s)/gamma(1 - a - s),
            (0, -re(a)/2 + S(1)/2), True)
    # TODO does not work with bneg, argument wrong. Needs changes to matching.
    #assert MT(expr.subs(b, -bpos), x, s) == \
    #       ((-1)**(a+1)*2**(a + 2*s)*bpos**(a + 2*s - 1)*gamma(a + s) \
    #                                      *gamma(1 - a - 2*s)/gamma(1 - s),
    #        (-re(a), -re(a)/2 + S(1)/2), True)
    expr = (sqrt(x+b**2)+b)**a
    assert MT(expr.subs(b, bpos), x, s) == \
           (-2**(a + 2*s)*a*bpos**(a + 2*s)*gamma(s)*gamma(-a - 2*s)/gamma(-a - s + 1),
            (0, -re(a)/2), True)
    #assert MT(expr.subs(b, -bpos), x, s) == \
    #       (2**(a + 2*s)*a*bpos**(a + 2*s)*gamma(-a - 2*s)*gamma(a + s)/gamma(-s + 1),
    #        (-re(a), -re(a)/2), True)
    # Test exponent 1:
    #assert MT(expr.subs({b: -bpos, a:1}), x, s) == \
    #       (-bpos**(2*s + 1)*gamma(s)*gamma(-s - S(1)/2)/(2*sqrt(pi)),
    #        (-1, -S(1)/2), True)

    # 8.4.2
    assert MT(exp(-x), x, s) == (gamma(s), (0, oo), True)
    assert MT(exp(-1/x), x, s) == (gamma(-s), (-oo, 0), True)

    # 8.4.5
    assert MT(log(x)**4*Heaviside(1-x), x, s) == (24/s**5, (0, oo), True)
    assert MT(log(x)**3*Heaviside(x-1), x, s) == (6/s**4, (-oo, 0), True)
    assert MT(log(x + 1), x, s) == (pi/(s*sin(pi*s)), (-1, 0), True)
    assert MT(log(1/x + 1), x, s) == (pi/(s*sin(pi*s)), (0, 1), True)
    assert MT(log(abs(1 - x)), x, s) == (pi/(s*tan(pi*s)), (-1, 0), True)
    assert MT(log(abs(1 - 1/x)), x, s) == (pi/(s*tan(pi*s)), (0, 1), True)

    # TODO we cannot currently do these (needs summation of 3F2(-1))
    #      this also implies that they cannot be written as a single g-function
    #      (although this is possible)
    mt = MT(log(x)/(x+1), x, s)
    assert mt[1:] == ((0, 1), True)
    assert not hyperexpand(mt[0], allow_hyper=True).has(meijerg)
    mt = MT(log(x)**2/(x+1), x, s)
    assert mt[1:] == ((0, 1), True)
    assert not hyperexpand(mt[0], allow_hyper=True).has(meijerg)
    mt = MT(log(x)/(x+1)**2, x, s)
    assert mt[1:] == ((0, 2), True)
    assert not hyperexpand(mt[0], allow_hyper=True).has(meijerg)

    # 8.4.14
    assert MT(erf(sqrt(x)), x, s) == \
           (-gamma(s + S(1)/2)/(sqrt(pi)*s), (-S(1)/2, 0), True)

    # 8.4.19
    assert MT(besselj(a, 2*sqrt(x)), x, s) == \
           (gamma(a/2 + s)/gamma(a/2 - s + 1), (-re(a)/2, S(3)/4), True)
    assert MT(sin(sqrt(x))*besselj(a, sqrt(x)), x, s) == \
           (2**a*gamma(S(1)/2 - 2*s)*gamma((a+1)/2 + s) \
                / (gamma(1 - s- a/2)*gamma(1 + a - 2*s)),
            (-(re(a) + 1)/2, S(1)/4), True)
    assert MT(cos(sqrt(x))*besselj(a, sqrt(x)), x, s) == \
           (2**a*gamma(a/2 + s)*gamma(S(1)/2 - 2*s)
                / (gamma(S(1)/2 - s - a/2)*gamma(a - 2*s + 1)),
            (-re(a)/2, S(1)/4), True)
    assert MT(besselj(a, sqrt(x))**2, x, s) == \
           (gamma(a + s)*gamma(S(1)/2 - s)
                / (sqrt(pi)*gamma(1 - s)*gamma(1 + a - s)),
            (-re(a), S(1)/2), True)
    assert MT(besselj(a, sqrt(x))*besselj(-a, sqrt(x)), x, s) == \
           (gamma(s)*gamma(S(1)/2 - s)
                / (sqrt(pi)*gamma(1 - a - s)*gamma(1 + a - s)),
            (0, S(1)/2), True)
    # NOTE: prudnikov gives the strip below as (1/2 - re(a), 1). As far as
    #       I can see this is wrong (since besselj(z) ~ 1/sqrt(z) for z large)
    assert MT(besselj(a - 1, sqrt(x))*besselj(a, sqrt(x)), x, s) == \
           (gamma(1-s)*gamma(a + s - S(1)/2)
                / (sqrt(pi)*gamma(S(3)/2 - s)*gamma(a - s + S(1)/2)),
            (S(1)/2 - re(a), S(1)/2), True)
    assert MT(besselj(a, sqrt(x))*besselj(b, sqrt(x)), x, s) == \
           (4**s*gamma(1 - 2*s)*gamma((a+b)/2 + s)
                / (gamma(1 - s + (b-a)/2)*gamma(1 - s + (a-b)/2)
                   *gamma( 1 - s + (a+b)/2)),
            (-(re(a) + re(b))/2, S(1)/2), True)
    assert MT(besselj(a, sqrt(x))**2 + besselj(-a, sqrt(x))**2, x, s)[1:] == \
           ((Max(re(a), -re(a)), S(1)/2), True)

    # Section 8.4.20
    assert MT(bessely(a, 2*sqrt(x)), x, s) == \
           (-cos(pi*a/2 - pi*s)*gamma(s - a/2)*gamma(s + a/2)/pi,
            (Max(-re(a)/2, re(a)/2), S(3)/4), True)
    assert MT(sin(sqrt(x))*bessely(a, sqrt(x)), x, s) == \
           (-4**s*sin(pi*a/2 - pi*s)*gamma(S(1)/2 - 2*s)
                * gamma((1-a)/2 + s)*gamma((1+a)/2 + s)
                / (sqrt(pi)*gamma(1 - s - a/2)*gamma(1 - s + a/2)),
            (Max(-(re(a) + 1)/2, (re(a) - 1)/2), S(1)/4), True)
    assert MT(cos(sqrt(x))*bessely(a, sqrt(x)), x, s) == \
           (-4**s*cos(pi*a/2 - pi*s)*gamma(s - a/2)*gamma(s + a/2)*gamma(S(1)/2 - 2*s)
                / (sqrt(pi)*gamma(S(1)/2 - s - a/2)*gamma(S(1)/2 - s + a/2)),
            (Max(-re(a)/2, re(a)/2), S(1)/4), True)
    assert MT(besselj(a, sqrt(x))*bessely(a, sqrt(x)), x, s) == \
           (-cos(pi*s)*gamma(s)*gamma(a + s)*gamma(S(1)/2 - s)
                / (pi**S('3/2')*gamma(1 + a - s)),
            (Max(-re(a), 0), S(1)/2), True)
    assert MT(besselj(a, sqrt(x))*bessely(b, sqrt(x)), x, s) == \
           (-4**s*cos(pi*a/2 - pi*b/2 + pi*s)*gamma(1 - 2*s)
                * gamma(a/2 - b/2 + s)*gamma(a/2 + b/2 + s)
                / (pi*gamma(a/2 - b/2 - s + 1)*gamma(a/2 + b/2 - s + 1)),
            (Max((-re(a) + re(b))/2, (-re(a) - re(b))/2), S(1)/2), True)
    assert MT(bessely(a, sqrt(x))**2, x, s) == \
           ((2*cos(pi*a - pi*s)*gamma(s)*gamma(-a + s)*gamma(-s + 1)*gamma(a - s + 1) \
             + pi*gamma(-s + S(1)/2)*gamma(s + S(1)/2)) \
            *gamma(a + s)/(pi**(S(3)/2)*gamma(-s + 1)*gamma(s + S(1)/2)*gamma(a - s + 1)),
            (Max(-re(a), 0, re(a)), S(1)/2), True)
    # TODO bessely(a, sqrt(x))*bessely(b, sqrt(x)) is a mess
    #      (no matter what way you look at it ...)

    # Section 8.4.22
    # TODO we can't do any of these (delicate cancellation)

    # Section 8.4.23
    assert MT(besselk(a, 2*sqrt(x)), x, s) == \
           (gamma(s - a/2)*gamma(s + a/2)/2, (Max(-re(a)/2, re(a)/2), oo), True)
    # TODO this result needs expansion of F(a, b; c; 1) using gauss-summation
    assert MT(exp(-x/2)*besselk(a, x/2), x, s)[1:] == \
           ((Max(-re(a), re(a)), oo), True)
    # TODO exp(x/2)*besselk(a, x/2) [etc] cannot currently be done
    # TODO products of bessel functions cannot be derived

def test_inverse_mellin_transform():
    from sympy import (sin, simplify, expand_func, powsimp, Max, Min, expand,
                       powdenest, powsimp, exp_polar)
    IMT = inverse_mellin_transform

    assert IMT(gamma(s), s, x, (0, oo)) == exp(-x)
    assert IMT(gamma(-s), s, x, (-oo, 0)) == exp(-1/x)
    assert IMT(s/(2*s**2 - 2), s, x, (2, oo)) \
           == (x**2/2 + S(1)/2)*Heaviside(1 - x)/(2*x)

    # test passing "None"
    assert IMT(1/(s**2 - 1), s, x, (-1, None)) \
           == -x*Heaviside(-x + 1)/2 - Heaviside(x - 1)/(2*x)
    assert IMT(1/(s**2 - 1), s, x, (None, 1)) \
           == -x*Heaviside(-x + 1)/2 - Heaviside(x - 1)/(2*x)

    # test expansion of sums
    assert IMT(gamma(s) + gamma(s-1), s, x, (1, oo)) == (x + 1)*exp(-x)/x

    # test factorisation of polys
    assert simplify(expand_func(IMT(1/(s**2 + 1), s, exp(-x),
                                    (None, oo))).rewrite(sin)) \
           == sin(x)*Heaviside(1 - exp(-x))

    # test multiplicative substitution
    a, b = symbols('a b', positive=True)
    c, d = symbols('c d')
    assert IMT(b**(-s/a)*factorial(s/a)/s, s, x, (0, oo)) == exp(-b*x**a)
    assert IMT(factorial(a/b + s/b)/(a+ s), s, x, (-a, oo)) == x**a*exp(-x**b)

    from sympy import expand_mul
    def simp_pows(expr): return expand_mul(simplify(powsimp(expr, force=True)), deep=True).replace(exp_polar, exp) # XXX ?

    # Now test the inverses of all direct transforms tested above

    # Section 8.4.2
    assert IMT(-1/(nu + s), s, x, (-oo, None)) == x**nu*Heaviside(x - 1)
    assert IMT(1/(nu + s), s, x, (None, oo)) == x**nu*Heaviside(1 - x)
    assert simp_pows(IMT(gamma(beta)*gamma(s)/gamma(s + beta), s, x, (0, oo))) \
           == (1 - x)**(beta - 1)*Heaviside(1 - x)
    assert simp_pows(IMT(gamma(beta)*gamma(1-beta-s)/gamma(1-s),
                         s, x, (-oo, None))) \
           == (x - 1)**(beta - 1)*Heaviside(x - 1)
    assert simp_pows(IMT(gamma(s)*gamma(rho-s)/gamma(rho), s, x, (0, None))) \
           == (1/(x + 1))**rho
    # TODO should this simplify further?
    assert simp_pows(IMT(d**c*d**(s-1)*sin(pi*c) \
                         *gamma(s)*gamma(s+c)*gamma(1-s)*gamma(1-s-c)/pi,
                         s, x, (Max(-re(c), 0), Min(1 - re(c), 1)))) \
           == d**c*(x/d)**c/(-d + x) - d**c/(-d + x)

    # TODO is calling simplify twice a bug?
    assert simplify(simplify(IMT(1/sqrt(pi)*(-c/2)*gamma(s)*gamma((1-c)/2 - s) \
                                 *gamma(-c/2-s)/gamma(1-c-s),
                                 s, x, (0, -re(c)/2)))) == \
           (1 + sqrt(x + 1))**c
    assert simplify(IMT(2**(a + 2*s)*b**(a + 2*s - 1)*gamma(s)*gamma(1 - a - 2*s) \
                        /gamma(1 - a - s), s, x, (0, (-re(a) + 1)/2))) == \
           (b + sqrt(b**2 + x))**(a - 1)*(b**2 + b*sqrt(b**2 + x) + x)/(b**2 + x)
    assert simplify(IMT(-2**(c + 2*s)*c*b**(c + 2*s)*gamma(s)*gamma(-c - 2*s) \
                          / gamma(-c - s + 1), s, x, (0, -re(c)/2))) == \
           (b + sqrt(b**2 + x))**c

    # Section 8.4.5
    assert IMT(24/s**5, s, x, (0, oo)) == log(x)**4*Heaviside(1 - x)
    assert expand(IMT(6/s**4, s, x, (-oo, 0)), force=True) == \
           log(x)**3*Heaviside(x - 1)
    assert IMT(pi/(s*sin(pi*s)), s, x, (-1, 0)) == log(x + 1)
    assert IMT(pi/(s*sin(pi*s/2)), s, x, (-2, 0)) == log(x**2 + 1)
    assert IMT(pi/(s*sin(2*pi*s)), s, x, (-S(1)/2, 0)) == log(sqrt(x) + 1)
    assert IMT(pi/(s*sin(pi*s)), s, x, (0, 1)) == log(1 + 1/x)

    # TODO
    def mysimp(expr):
        from sympy import expand, logcombine, powsimp
        return expand(powsimp(logcombine(expr, force=True), force=True, deep=True),
                      force=True).replace(exp_polar, exp)
    assert mysimp(mysimp(IMT(pi/(s*tan(pi*s)), s, x, (-1, 0)))) == \
           log(1-x)*Heaviside(1-x) + log(x-1)*Heaviside(x-1)
    assert mysimp(IMT(pi/(s*tan(pi*s)), s, x, (0, 1))) == \
           log(1/x - 1)*Heaviside(1-x) + log(1 - 1/x)*Heaviside(x-1)

    # 8.4.14
    assert IMT(-gamma(s + S(1)/2)/(sqrt(pi)*s), s, x, (-S(1)/2, 0)) == \
           erf(sqrt(x))

    # 8.4.19
    # TODO these come out ugly
    def mysimp(expr):
        return expand(unpolarify(simplify(expand(expand_func(expr.rewrite(besselj))))))
    assert mysimp(IMT(gamma(a/2 + s)/gamma(a/2 - s + 1), s, x, (-re(a)/2, S(3)/4))) \
           == besselj(a, 2*sqrt(x)*polar_lift(-1))*exp(-I*pi*a)
    assert mysimp(IMT(2**a*gamma(S(1)/2 - 2*s)*gamma(s + (a + 1)/2) \
                      / (gamma(1 - s - a/2)*gamma(1 - 2*s + a)),
                      s, x, (-(re(a) + 1)/2, S(1)/4))) == \
           exp(-I*pi*a)*sin(sqrt(x))*besselj(a, sqrt(x)*polar_lift(-1))
    assert mysimp(IMT(2**a*gamma(a/2 + s)*gamma(S(1)/2 - 2*s) \
                      / (gamma(S(1)/2 - s - a/2)*gamma(1 - 2*s + a)),
                      s, x, (-re(a)/2, S(1)/4))) == \
           exp(-I*pi*a)*cos(sqrt(x))*besselj(a, sqrt(x)*polar_lift(-1))
    # TODO this comes out as an amazing mess, but surprisingly enough mysimp is
    #      effective ...
    assert powsimp(powdenest(mysimp(IMT(gamma(a + s)*gamma(S(1)/2 - s) \
                      / (sqrt(pi)*gamma(1 - s)*gamma(1 + a - s)),
                      s, x, (-re(a), S(1)/2))), polar=True)) == \
           exp(-2*I*pi*a)*besselj(a, sqrt(x)*polar_lift(-1))**2
    # NOTE the next is indeed an even function of sqrt(x), so the result is
    #      correct
    assert mysimp(IMT(gamma(s)*gamma(S(1)/2 - s) \
                      / (sqrt(pi)*gamma(1 - s - a)*gamma(1 + a - s)),
                      s, x, (0, S(1)/2))) == \
           besselj(-a, polar_lift(-1)*sqrt(x))*besselj(a, polar_lift(-1)*sqrt(x))
    assert mysimp(IMT(4**s*gamma(-2*s + 1)*gamma(a/2 + b/2 + s) \
                      / (gamma(-a/2 + b/2 - s + 1)*gamma(a/2 - b/2 - s + 1) \
                         *gamma(a/2 + b/2 - s + 1)),
                      s, x, (-(re(a) + re(b))/2, S(1)/2))) == \
           exp(-I*pi*a)*exp(-I*pi*b)*besselj(a, sqrt(x)*polar_lift(-1)) \
                                    *besselj(b, sqrt(x)*polar_lift(-1))

    # Section 8.4.20
    # TODO these come out even messier, not worth testing for now

    # TODO the other bessel functions, when simplification is there

def test_laplace_transform():
    LT = laplace_transform
    a, b, c, = symbols('a b c', positive=True)
    t = symbols('t')

    # basic tests from wikipedia

    assert LT((t-a)**b*exp(-c*(t-a))*Heaviside(t-a), t, s) \
           == ((s + c)**(-b - 1)*exp(-a*s)*gamma(b + 1), -c, True)
    assert LT(t**a, t, s) == (s**(-a - 1)*gamma(a + 1), 0, True)
    assert LT(Heaviside(t), t, s) == (1/s, 0, True)
    assert LT(Heaviside(t - a), t, s) == (exp(-a*s)/s, 0, True)
    assert LT(1 - exp(-a*t), t, s) == (a/(s*(a + s)), 0, True)

    assert LT((exp(2*t)-1)*exp(-b - t)*Heaviside(t)/2, t, s, noconds=True) \
           == exp(-b)/(s**2 - 1)

    assert LT(exp(t), t, s)[0:2] == (1/(s-1), 1)
    assert LT(exp(2*t), t, s)[0:2] == (1/(s-2), 2)
    assert LT(exp(a*t), t, s)[0:2] == (1/(s-a), a)

    lt = LT(log(t/a), t, s)
    assert lt[1:] == (0, True)
    # TODO hyperexpand is not clever enough to recognise this on its own
    assert hyperexpand(lt[0], allow_hyper=True) == (-log(a*s) - EulerGamma)/s

    assert LT(erf(t), t, s) == ((-erf(s/2) + 1)*exp(s**2/4)/s, 0, True)

    assert LT(sin(a*t), t, s) == (a/(a**2 + s**2), 0, True)
    assert LT(cos(a*t), t, s) == (s/(a**2 + s**2), 0, True)
    # TODO would be nice to have these come out better
    assert LT(exp(-a*t)*sin(b*t), t, s) == (1/b/(1 + (a + s)**2/b**2), -a, True)
    assert LT(exp(-a*t)*cos(b*t), t, s) == \
           (sqrt((a + s)**2)/b**2/(1 + (a + s)**2/b**2), -a, True)
    # TODO sinh, cosh have delicate cancellation

    # TODO conditions are a mess
    assert LT(besselj(0, t), t, s, noconds=True) == 1/sqrt(1 + s**2)
    assert LT(besselj(1, t), t, s, noconds=True) == 1 - 1/sqrt(1 + 1/s**2)
    # TODO general order works, but is a *mess*
    # TODO besseli also works, but is an eaven greater mess

def test_inverse_laplace_transform():
    from sympy import expand, sinh, cosh, besselj, besseli, exp_polar, unpolarify
    ILT = inverse_laplace_transform
    a, b, c, = symbols('a b c', positive=True)
    t = symbols('t')

    def simp_hyp(expr):
        return expand(expand(expr).rewrite(sin))

    # just test inverses of all of the above
    assert ILT(1/s, s, t) == Heaviside(t)
    assert ILT(1/s**2, s, t) == t*Heaviside(t)
    assert ILT(1/s**5, s, t) == t**4*Heaviside(t)/factorial(4)
    assert ILT(exp(-a*s)/s, s, t) == Heaviside(t-a)
    assert ILT(exp(-a*s)/(s+b), s, t) == exp(b*(a - t))*Heaviside(t - a)
    assert ILT(a/(s**2 + a**2), s, t) == sin(a*t)*Heaviside(t)
    assert ILT(s/(s**2 + a**2), s, t) == cos(a*t)*Heaviside(t)
    # TODO is there a way around simp_hyp?
    assert simp_hyp(ILT(a/(s**2 - a**2), s, t)) == sinh(a*t)*Heaviside(t)
    assert simp_hyp(ILT(s/(s**2 - a**2), s, t)) == cosh(a*t)*Heaviside(t)
    assert ILT(a/((s+b)**2 + a**2), s, t) == exp(-b*t)*sin(a*t)*Heaviside(t)
    assert ILT((s+b)/((s+b)**2 + a**2), s, t) == exp(-b*t)*cos(a*t)*Heaviside(t)
    # TODO sinh/cosh shifted come out a mess. also delayed trig is a mess
    # TODO should this simplify further?
    assert ILT(exp(-a*s)/s**b, s, t) == \
      (t - a)**(b - 1)*Heaviside(t - a)/gamma(b)

    assert ILT(exp(-a*s)/sqrt(1 + s**2), s, t) == \
        Heaviside(t - a)*besselj(0, a - t) # note: besselj(0, x) is even

    # TODO besselsimp would be good to have
    # these are re-enabled in a later commit
    #assert ILT(a**b*(s + sqrt(s**2 - a**2))**(-b)/sqrt(s**2 - a**2), s, t) == \
    #    (t**2)**(b/2 + S(1)/2)*Heaviside(t)*besseli(b, a*t*exp_polar(I*pi))*exp_polar(-I*pi*b)/(t*sqrt((t**2)**b))
    #assert unpolarify(ILT(a**b*(s + sqrt(s**2 + a**2))**(-b)/sqrt(s**2 + a**2),
    #                      s, t).rewrite(besselj)) == \
    #    (t**2)**(b/2 + S(1)/2)*exp(-I*pi*b)*Heaviside(t)*besselj(b, a*t*exp_polar(I*pi))/(t*sqrt((t**2)**b))

    assert ILT(1/(s*sqrt(s+1)), s, t) == Heaviside(t)*erf(sqrt(t))
    # TODO can we make erf(t) work?

def test_fourier_transform():
    from sympy import simplify, expand, expand_complex, factor, expand_trig
    FT = fourier_transform
    IFT = inverse_fourier_transform
    def simp(x): return simplify(expand_trig(expand_complex(expand(x))))
    def sinc(x): return sin(pi*x)/(pi*x)
    k = symbols('k', real=True)

    # TODO for this to work with real a, need to expand abs(a*x) to abs(a)*abs(x)
    a = symbols('a', positive=True)
    b = symbols('b', positive=True)

    posk = symbols('k', positive=True)

    # basic examples from wikipedia
    assert simp(FT(Heaviside(1 - abs(2*a*x)), x, k)) == sinc(k/a)/a
    # TODO IFT is a *mess*
    assert simp(FT(Heaviside(1-abs(a*x))*(1-abs(a*x)), x, k)) == sinc(k/a)**2/a
    # TODO IFT

    assert factor(FT(exp(-a*x)*Heaviside(x), x, k), extension=I) \
           == 1/(a + 2*pi*I*k)
    # NOTE: the ift comes out in pieces
    assert IFT(1/(a + 2*pi*I*x), x, posk, noconds=False) == (exp(-a*posk), True)
    assert IFT(1/(a + 2*pi*I*x), x, -posk,
               noconds=False) == (0, True)
    assert IFT(1/(a + 2*pi*I*x), x, symbols('k', negative=True),
               noconds=False) == (0, True)
    # TODO IFT without factoring comes out as meijer g

    assert factor(FT(x*exp(-a*x)*Heaviside(x), x, k), extension=I) \
           == 1/(a + 2*pi*I*k)**2
    assert FT(exp(-a*x)*sin(b*x)*Heaviside(x), x, k) \
           == 1/b/(1 + (a + 2*I*pi*k)**2/b**2)

    assert FT(exp(-a*x**2), x, k) == sqrt(pi)*exp(-pi**2*k**2/a)/sqrt(a)
    assert IFT(sqrt(pi/a)*exp(-(pi*k)**2/a), k, x) == exp(-a*x**2)
    assert FT(exp(-a*abs(x)), x, k) == 2*a/(a**2 + 4*pi**2*k**2)
    # TODO IFT (comes out as meijer G)

    # TODO besselj(n, x), n an integer > 0 actually can be done...

    # TODO are there other common transforms (no distributions!)?
