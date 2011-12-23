from sympy import (meijerg, I, S, integrate, Integral, oo, gamma,
                   hyperexpand, exp, simplify, sqrt, pi, erf, sin, cos)
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

    # Test trigexpand:
    assert integrate(exp(-x)*sin(x + a), (x, 0, oo), meijerg=True) == \
           sin(a)/2 + cos(a)/2

def test_bessel():
    from sympy import besselj, Heaviside, besseli, polar_lift, exp_polar
    assert integrate(besselj(a, z)*besselj(b, z)/z, (z, 0, oo),
                     meijerg=True, conds='none') == \
           2*sin(pi*a/2 - pi*b/2)/(pi*(a-b)*(a+b))
    assert integrate(besselj(a, z)*besselj(a, z)/z, (z, 0, oo),
                     meijerg=True, conds='none') == 1/(2*a)

    # TODO more orthogonality integrals

    # TODO there is actually a lot to improve here, this example is a good
    #      stress-test
    # (the original integral is not recognised, the convergence conditions are
    #  wrong, and the result can be simplified to besselj(y, z))
    assert simplify(integrate(sin(z*x)*(x**2-1)**(-(y+S(1)/2))*Heaviside(x**2-1),
                              (x, 0, oo), meijerg=True, conds='none')
                    *2/((z/2)**y*sqrt(pi)*gamma(S(1)/2-y))) == \
           2*(z**2/4)**(y + S(1)/2)*(z/2)**(-y)*(2*exp_polar(-I*pi/2)/z)**y \
           *besseli(y, z*exp_polar(I*pi/2))/z

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
        for formula, terms, cond in sorted(l):
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
