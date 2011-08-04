from sympy import (meijerg, I, S, integrate, Integral, oo, gamma,
                   hyperexpand, exp, simplify, sqrt, pi, erf, sin, cos,
                   exp_polar, polar_lift, polygamma, hyper, log)
from sympy.integrals.meijerint import (_rewrite_single, _rewrite1,
         meijerint_indefinite, _inflate_g, _create_lookup_table,
         meijerint_definite, meijerint_inversion)
from sympy.utilities.randtest import (test_numerically,
         random_complex_number as randcplx)
from sympy.abc import x, y, a, b, c, d, s, t, z

def test_rewrite_single():
    def t(expr, c, m):
        e = _rewrite_single(meijerg([a], [b], [c], [d], expr), x)
        assert e is not None
        assert isinstance(e[0][0][2], meijerg)
        assert e[0][0][2].argument.as_coeff_mul(x) == (c, (m,))
    def tn(expr):
        assert _rewrite_single(meijerg([a], [b], [c], [d], expr), x) is None

    t(x, 1, x)
    t(x**2, 1, x**2)
    t(x**2 + y*x**2, y + 1, x**2)
    tn(x**2 + x)
    tn(x**y)

    def u(expr, x):
        from sympy import Add, exp, exp_polar
        r = _rewrite_single(expr, x)
        e = Add(*[res[0]*res[2] for res in r[0]]).replace(exp_polar, exp) # XXX Hack?
        assert test_numerically(e, expr, x)

    u(exp(-x)*sin(x), x)
    u(exp(-x)*sin(x)*cos(x), x)
    u(exp(x)*sin(x), x)

def test_rewrite1():
    assert _rewrite1(x**3*meijerg([a], [b], [c], [d], x**2 + y*x**2)*5, x) \
           == (5, x**3, [(1, 0, meijerg([a], [b], [c], [d], x**2*(y + 1)))], True)

def test_meijerint_indefinite_numerically():
    def t(fac, arg):
        g = meijerg([a], [b], [c], [d], arg)*fac
        subs = {a: randcplx()/10, b:randcplx()/10 + I,
                c: randcplx(), d: randcplx()}
        integral = meijerint_indefinite(g, x)
        assert integral is not None
        assert test_numerically(g.subs(subs), integral.diff(x).subs(subs), x)
    t(1, x)
    t(2, x)
    t(1, 2*x)
    t(1, x**2)
    t(5, x**S('3/2'))
    t(x**3, x)
    t(3*x**S('3/2'), 4*x**S('7/3'))

def test_inflate():
    subs = {a: randcplx()/10, b: randcplx()/10 + I, c: randcplx(),
            d: randcplx(), y:randcplx()/10}
    def t(a, b, arg, n):
        from sympy import Mul
        m1 = meijerg(a, b, arg)
        m2 = Mul(*_inflate_g(m1, n))
        # NOTE: (the random number)**9 must still be on the principal sheet.
        # Thus make b, d small to create random numbers of small imaginary part.
        return test_numerically(m1.subs(subs), m2.subs(subs), x, b=0.1, d=-0.1)
    assert t([[a], [b]], [[c], [d]], x, 3)
    assert t([[a, y], [b]], [[c], [d]], x, 3)
    assert t([[a], [b]], [[c, y], [d]], 2*x**3, 3)

def test_recursive():
    from sympy import symbols, exp_polar
    a, b, c = symbols('a b c', positive=True)
    assert simplify(integrate(exp(-(x-a)**2)*exp(-(x-b)**2), (x, 0, oo))) \
           == sqrt(2*pi)/4*(1 + erf(sqrt(2)/2*(a + b))) \
              *exp(-a**2 - b**2 + (a + b)**2/2)
    assert simplify(integrate(exp(-(x-a)**2)*exp(-(x-b)**2)*exp(c*x), (x, 0, oo))) \
           == sqrt(2*pi)/4*(1 + erf(sqrt(2)/4*(2*a + 2*b + c))) \
              *exp(-a**2 - b**2 + (2*a + 2*b + c)**2/8)
    assert simplify(integrate(exp(-(x-a-b-c)**2), (x, 0, oo))) \
           == sqrt(pi)/2*(1 + erf(a + b + c))
    assert simplify(integrate(exp(-(x+a+b+c)**2), (x, 0, oo))) \
           == sqrt(pi)/2*(1 - erf(a + b + c))

def test_meijerint():
    from sympy import symbols, expand, arg
    s, t, mu = symbols('s t mu', real=True)
    assert integrate(meijerg([], [], [0], [], s*t)
                     *meijerg([], [], [mu/2], [-mu/2], t**2/4),
                     (t, 0, oo)).is_Piecewise
    s = symbols('s', positive=True)
    assert integrate(x**s*meijerg([[],[]], [[0],[]], x), (x, 0, oo)) \
           == gamma(s + 1)
    assert integrate(x**s*meijerg([[],[]], [[0],[]], x), (x, 0, oo),
                     meijerg=True) == gamma(s + 1)
    assert isinstance(integrate(x**s*meijerg([[],[]], [[0],[]], x),
                                (x, 0, oo), meijerg=False),
                      Integral)

    assert meijerint_indefinite(exp(x), x) == exp(x)

    # TODO what simplifications should be done automatically?
    # This tests "extra case" for antecedents_1.
    a, b = symbols('a b', positive=True)
    assert simplify(meijerint_definite(x**a, x, 0, b)[0]) \
           == b**(a + 1)/(a + 1)

    # This tests various conditions and expansions:
    meijerint_definite((x+1)**3*exp(-x), x, 0, oo) == (16, True)

    # Again, how about simplifications?
    sigma, mu = symbols('sigma mu', positive=True)
    i, c = meijerint_definite(exp(-((x-mu)/(2*sigma))**2), x, 0, oo)
    assert simplify(i) \
           == sqrt(pi)*sigma*(erf(mu/(2*sigma)) + 1)
    assert c is True

    i, _ = meijerint_definite(exp(-mu*x)*exp(sigma*x), x, 0, oo)
    # TODO it would be nice to test the condition
    assert simplify(i) == 1/(mu - sigma)

    # Test substitutions to change limits
    assert meijerint_definite(exp(x), x, -oo, 2) == (exp(2), True)
    assert expand(meijerint_definite(exp(x), x, 0, I)[0]) == exp(I) - 1
    assert expand(meijerint_definite(exp(-x), x, 0, x)[0]) == \
           1 - exp(-exp(I*arg(x))*abs(x))

    # Test -oo to oo
    assert meijerint_definite(exp(-x**2), x, -oo, oo) == (sqrt(pi), True)
    assert meijerint_definite(exp(-abs(x)), x, -oo, oo) == (2, True)
    assert meijerint_definite(exp(-(2*x-3)**2), x, -oo, oo) == (sqrt(pi)/2, True)
    assert meijerint_definite(exp(-abs(2*x-3)), x, -oo, oo) == (1, True)
    assert meijerint_definite(exp(-((x-mu)/sigma)**2/2)/sqrt(2*pi*sigma**2),
                              x, -oo, oo) == (1, True)

    # Test one of the extra conditions for 2 g-functinos
    assert meijerint_definite(exp(-x)*sin(x), x, 0, oo) == (S(1)/2, True)

    # Test a bug
    def res(n): return (1/(1+x**2)).diff(x, n).subs(x,1)*(-1)**n
    for n in range(6):
       assert integrate(exp(-x)*sin(x)*x**n, (x, 0, oo), meijerg=True) == res(n)

    # This used to test trigexpand... now it is done by linear substitution
    assert simplify(integrate(exp(-x)*sin(x + a), (x, 0, oo), meijerg=True)).expand().rewrite(sin).expand() == \
           sin(a)/2 + cos(a)/2

    # This is cosine-integral, not currently in hyperexpand tables
    from sympy import EulerGamma
    a = symbols('a', positive=True)
    assert simplify(integrate(cos(x)/x, (x, a, oo), meijerg=True)) == \
           a**2*hyper((1, 1), (S(3)/2, 2, 2), a**2*exp_polar(I*pi)/4)/4 \
           + log(2/a) + polygamma(0, S(1)/2)/2 - EulerGamma/2

    # Test the condition 14 from prudnikov.
    # (This is besselj*besselj in disguise, to stop the product from being
    #  recognised in the tables.)
    a, b, s = symbols('a b s')
    from sympy import And, re
    assert meijerint_definite(meijerg([], [], [a/2], [-a/2], x/4) \
                  *meijerg([], [], [b/2], [-b/2], x/4)*x**(s-1), x, 0, oo) == \
           (4*2**(2*s - 2)*gamma(-2*s + 1)*gamma(a/2 + b/2 + s) \
               /(gamma(-a/2 + b/2 - s + 1)*gamma(a/2 - b/2 - s + 1) \
                 *gamma(a/2 + b/2 - s + 1)),
            And(0 < -2*re(4*s) + 8, 0 < re(a/2 + b/2 + s), re(2*s) < 1))

    # test a bug
    assert integrate(sin(x**a)*sin(x**b), (x, 0, oo), meijerg=True) == \
           Integral(sin(x**a)*sin(x**b), (x, 0, oo))

    # test better hyperexpand
    assert integrate(exp(-x**2)*log(x), (x, 0, oo), meijerg=True) == \
           sqrt(pi)*polygamma(0, S(1)/2)/4
    n = symbols('n', integer = True)
    assert simplify(integrate(exp(-x)*x**n, x, meijerg=True)) == \
           x**(n + 1)*hyper([n + 1], [n + 2], polar_lift(-1)*x)/(n + 1)

def test_bessel():
    from sympy import (besselj, Heaviside, besseli, polar_lift, exp_polar,
                       powdenest)
    assert simplify(integrate(besselj(a, z)*besselj(b, z)/z, (z, 0, oo),
                     meijerg=True, conds='none')) == \
           2*sin(pi*a/2 - pi*b/2)/(pi*(a-b)*(a+b))
    assert simplify(integrate(besselj(a, z)*besselj(a, z)/z, (z, 0, oo),
                     meijerg=True, conds='none')) == 1/(2*a)

    # TODO more orthogonality integrals

    # TODO there is some improvement possible here:
    #  - the result can be simplified to besselj(y, z))
    assert powdenest(simplify(integrate(sin(z*x)*(x**2-1)**(-(y+S(1)/2)),
                              (x, 1, oo), meijerg=True, conds='none')
                              *2/((z/2)**y*sqrt(pi)*gamma(S(1)/2-y))),
                     polar=True) == \
           exp(-I*pi*y/2)*besseli(y, z*exp_polar(I*pi/2))

    # Werner Rosenheinrich
    # SOME INDEFINITE INTEGRALS OF BESSEL FUNCTIONS

    assert integrate(x*besselj(0, x), x, meijerg=True) == x*besselj(1, x)
    assert integrate(x*besseli(0, x), x, meijerg=True) == x*besseli(1, x)
    # TODO can do higher powers, but come out as high order ... should they be
    #      reduced to order 0, 1?
    assert integrate(besselj(1, x), x, meijerg=True) == -besselj(0, x)
    assert integrate(besselj(1, x)**2/x, x, meijerg=True) == \
           -(besselj(0, x)**2 + besselj(1, x)**2)/2
    # TODO more besseli when tables are extended or recursive mellin works
    assert integrate(besselj(0, x)**2/x**2, x, meijerg=True) == \
           -2*x*besselj(0, x)**2 - 2*x*besselj(1, x)**2 \
           + 2*besselj(0, x)*besselj(1, x) - besselj(0, x)**2/x
    assert integrate(besselj(0, x)*besselj(1, x), x, meijerg=True) == \
           -besselj(0, x)**2/2
    assert integrate(x**2*besselj(0, x)*besselj(1, x), x, meijerg=True) == \
           x**2*besselj(1, x)**2/2
    assert integrate(besselj(0, x)*besselj(1, x)/x, x, meijerg=True) == \
           x*besselj(0, x)**2 + x*besselj(1, x)**2 - besselj(0, x)*besselj(1, x)
    # TODO how does besselj(0, a*x)*besselj(0, b*x) work?
    # TODO how does besselj(0, x)**2*besselj(1, x)**2 work?
    # TODO sin(x)*besselj(0, x) etc come out a mess
    # TODO can x*log(x)*besselj(0, x) be done?
    # TODO how does besselj(1, x)*besselj(0, x+a) work?
    # TODO more indefinite integrals when struve functions etc are implemented

    # test a substitution
    assert integrate(besselj(1, x**2)*x, x, meijerg=True) == -besselj(0, x**2)/2

def test_inversion():
    from sympy import piecewise_fold, besselj, sqrt, I, sin, cos, Heaviside
    def inv(f): return piecewise_fold(meijerint_inversion(f, s, t))
    assert inv(1/(s**2 + 1)) == sin(t)*Heaviside(t)
    assert inv(s/(s**2 + 1)) == cos(t)*Heaviside(t)
    assert inv(exp(-s)/s) == Heaviside(t - 1)
    assert inv(1/sqrt(1 + s**2)) == besselj(0, t)*Heaviside(t)

    # Test some antcedents checking.
    assert meijerint_inversion(sqrt(s)/sqrt(1 + s**2), s, t) is None
    assert inv(exp(s**2)) is None
    assert meijerint_inversion(exp(-s**2), s, t) is None

def test_lookup_table():
    from random import uniform, randrange
    from sympy import Add, unpolarify, exp_polar, exp
    from sympy.integrals.meijerint import z as z_dummy
    table = {}
    _create_lookup_table(table)
    for _, l in sorted(table.items()):
        for formula, terms, cond, hint in sorted(l):
            subs = {}
            for a in list(formula.free_symbols) + [z_dummy]:
                if hasattr(a, 'properties') and a.properties:
                    # these Wilds match positive integers
                    subs[a] = randrange(1, 10)
                else:
                    subs[a] = uniform(1.5, 3.5)
            if not isinstance(terms, list):
                terms = terms(subs)

            # First test that hyperexpand can do this.
            expanded = [hyperexpand(g) for (_, g) in terms]
            assert all (x.is_Piecewise or not x.has(meijerg) for x in expanded)

            # Now test that the meijer g-function is indeed as advertised.
            expanded = Add(*[f*x for (f, x) in terms])
            a, b = formula.n(subs=subs), expanded.n(subs=subs)
            r = min(abs(a), abs(b))
            if r < 1:
                assert abs(a - b).n() <= 1e-10
            else:
                assert (abs(a - b)/r).n() <= 1e-10

def test_branch_bug():
    from sympy import powdenest, lowergamma
    # TODO combsimp cannot prove that the factor is unity
    assert powdenest(integrate(erf(x**3), x, meijerg=True).diff(x), polar=True) \
           == 2*erf(x**3)*gamma(S(2)/3)/3/gamma(S(5)/3)
    assert integrate(erf(x**3), x, meijerg=True) == \
           2*x*erf(x**3)*gamma(S(2)/3)/(3*gamma(S(5)/3)) \
           - 2*gamma(S(2)/3)*lowergamma(S(2)/3, x**6)/(3*sqrt(pi)*gamma(S(5)/3))

def test_linear_subs():
    from sympy import besselj
    assert integrate(sin(x-1), x, meijerg=True) == -cos(1 - x)
    assert integrate(besselj(1, x-1), x, meijerg=True) == -besselj(0, 1 - x)

def test_probability():
    # various integrals from probability theory
    from sympy.abc import x, y, z
    from sympy import symbols, Symbol, Abs, expand_mul, combsimp, powsimp
    mu1, mu2 = symbols('mu1 mu2', real=True, finite=True, bounded=True)
    sigma1, sigma2 = symbols('sigma1 sigma2', real=True, finite=True,
                                              bounded=True, positive=True)
    rate = Symbol('lambda', real=True, positive=True, bounded=True)
    def normal(x, mu, sigma):
        return 1/sqrt(2*pi*sigma**2)*exp(-(x-mu)**2/2/sigma**2)
    def exponential(x, rate):
        return rate*exp(-rate*x)

    assert integrate(normal(x, mu1, sigma1), (x, -oo, oo), meijerg=True) == 1
    assert integrate(x*normal(x, mu1, sigma1), (x, -oo, oo), meijerg=True) == mu1
    assert integrate(x**2*normal(x, mu1, sigma1), (x, -oo, oo), meijerg=True) \
           == mu1**2 + sigma1**2
    assert integrate(x**3*normal(x, mu1, sigma1), (x, -oo, oo), meijerg=True) \
           == mu1**3 + 3*mu1*sigma1**2
    assert integrate(normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                     (x, -oo, oo), (y, -oo, oo), meijerg=True) == 1
    assert integrate(x*normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                     (x, -oo, oo), (y, -oo, oo), meijerg=True) == mu1
    assert integrate(y*normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                     (x, -oo, oo), (y, -oo, oo), meijerg=True) == mu2
    assert integrate(x*y*normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                     (x, -oo, oo), (y, -oo, oo), meijerg=True) == mu1*mu2
    assert integrate((x+y+1)*normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                     (x, -oo, oo), (y, -oo, oo), meijerg=True) == 1 + mu1 + mu2
    assert integrate((x+y-1)*normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                     (x, -oo, oo), (y, -oo, oo), meijerg=True) == -1 + mu1 + mu2

    i = integrate(x**2*normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                  (x, -oo, oo), (y, -oo, oo), meijerg=True)
    assert not i.has(Abs)
    assert simplify(i) == mu1**2 + sigma1**2
    assert integrate(y**2*normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                     (x, -oo, oo), (y, -oo, oo), meijerg=True) == \
           sigma2**2 + mu2**2

    assert integrate(exponential(x, rate), (x, 0, oo), meijerg=True) == 1
    assert integrate(x*exponential(x, rate), (x, 0, oo), meijerg=True) == 1/rate
    assert integrate(x**2*exponential(x, rate), (x, 0, oo), meijerg=True) \
           == 2/rate**2

    def E(expr):
        res1 = integrate(expr*exponential(x, rate)*normal(y, mu1, sigma1),
                         (x, 0, oo), (y, -oo, oo), meijerg=True)
        res2 = integrate(expr*exponential(x, rate)*normal(y, mu1, sigma1),
                                 (y, -oo, oo), (x, 0, oo), meijerg=True)
        assert expand_mul(res1) == expand_mul(res2)
        return res1

    assert E(1) == 1
    assert E(x*y) == mu1/rate
    assert E(x*y**2) == mu1**2/rate + sigma1**2/rate
    assert simplify(E((x+y+1)**2) - E(x+y+1)**2) == (rate**2*sigma1**2 + 1)/rate**2
    assert simplify(E((x+y-1)**2) - E(x+y-1)**2) == (rate**2*sigma1**2 + 1)/rate**2
    assert simplify(E((x+y)**2) - E(x+y)**2) == (rate**2*sigma1**2 + 1)/rate**2

    # Beta' distribution
    alpha, beta = symbols('alpha beta', positive=True)
    betadist = x**(alpha-1)*(1+x)**(-alpha - beta)*gamma(alpha+beta) \
              /gamma(alpha)/gamma(beta)
    assert integrate(betadist, (x, 0, oo), meijerg=True) == 1
    i = integrate(x*betadist, (x, 0, oo), meijerg=True, conds='separate')
    assert (combsimp(i[0]), i[1]) == (alpha/(beta - 1), 1 < beta)
    j = integrate(x**2*betadist, (x, 0, oo), meijerg=True, conds='separate')
    assert j[1] == (1 < beta - 1)
    assert combsimp(j[0] - i[0]**2) == (alpha + beta - 1)*alpha \
                                        /(beta - 2)/(beta - 1)**2

    # Chi distribution
    k = Symbol('k', integer=True, positive=True)
    chi = 2**(1-k/2)*x**(k-1)*exp(-x**2/2)/gamma(k/2)
    assert powsimp(integrate(chi, (x, 0, oo), meijerg=True)) == 1
    assert simplify(integrate(x*chi, (x, 0, oo), meijerg=True)) == \
           sqrt(2)*gamma((k + 1)/2)/gamma(k/2)
    assert simplify(integrate(x**2*chi, (x, 0, oo), meijerg=True)) == k

    # Chi^2 distribution
    chisquared = 2**(-k/2)/gamma(k/2)*x**(k/2-1)*exp(-x/2)
    assert powsimp(integrate(chisquared, (x, 0, oo), meijerg=True)) == 1
    assert simplify(integrate(x*chisquared, (x, 0, oo), meijerg=True)) == k
    assert simplify(integrate(x**2*chisquared, (x, 0, oo), meijerg=True)) == \
           k*(k + 2)
    assert combsimp(integrate(((x-k)/sqrt(2*k))**3*chisquared, (x, 0, oo),
                    meijerg=True)) == 2*sqrt(2)/sqrt(k)

    # Dagum distribution
    a, b, p = symbols('a b p', positive=True)
    # XXX (x/b)**a does not work
    dagum = a*p/x*(x/b)**(a*p)/(1 + x**a/b**a)**(p+1)
    assert simplify(integrate(dagum, (x, 0, oo), meijerg=True)) == 1
    # XXX conditions are a mess
    assert simplify(integrate(x*dagum, (x, 0, oo), meijerg=True, conds='none')) \
           == b*gamma(1 - 1/a)*gamma(p + 1/a)/gamma(p)
    assert simplify(integrate(x**2*dagum, (x, 0, oo), meijerg=True, conds='none')) \
           == b**2*gamma(1 - 2/a)*gamma(p + 2/a)/gamma(p)

    # F-distribution
    d1, d2 = symbols('d1 d2', positive=True)
    f = sqrt(((d1*x)**d1 * d2**d2)/(d1*x + d2)**(d1+d2))/x \
          /gamma(d1/2)/gamma(d2/2)*gamma((d1 + d2)/2)
    assert simplify(integrate(f, (x, 0, oo), meijerg=True)) == 1
    # TODO conditions are a mess
    assert simplify(integrate(x*f, (x, 0, oo), meijerg=True, conds='none')) == \
           d2/(d2 - 2)
    assert simplify(integrate(x**2*f, (x, 0, oo), meijerg=True, conds='none')) == \
           d2**2*(d1 + 2)/d1/(d2 - 4)/(d2 - 2)

    # TODO gamma, inverse gaussian, Levi, log-logistic, rayleigh, weibull

    # rice distribution
    from sympy import besseli
    nu, sigma = symbols('nu sigma', positive=True)
    rice = x/sigma**2*exp(-(x**2+ nu**2)/2/sigma**2)*besseli(0, x*nu/sigma**2)
    assert integrate(rice, (x, 0, oo), meijerg=True) == 1
    # can someone verify higher moments?

    # Laplace distribution
    mu = Symbol('mu', real=True)
    b = Symbol('b', positive=True)
    laplace = exp(-abs(x-mu)/b)/2/b
    assert integrate(laplace, (x, -oo, oo), meijerg=True) == 1
    assert integrate(x*laplace, (x, -oo, oo), meijerg=True) == mu
    assert integrate(x**2*laplace, (x, -oo, oo), meijerg=True) == 2*b**2 + mu**2

    # TODO are there other distributions supported on (-oo, oo) that we can do?

    # misc tests
    k = Symbol('k', integer=True)
    assert simplify(integrate(log(x) * x**(k-1) * exp(-x) / gamma(k), (x, 0, oo),
                              conds='none')) == polygamma(0, k)

