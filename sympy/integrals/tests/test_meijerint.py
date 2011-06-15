from sympy import (meijerg, I, S, integrate, Integral, oo, gamma,
                   hyperexpand, exp, simplify, sqrt, pi, erf)
from sympy.integrals.meijerint import (_rewrite_single, _rewrite1,
         meijerint_indefinite, _inflate_g, _create_lookup_table,
         meijerint_definite)
from sympy.utilities.randtest import (test_numerically,
         random_complex_number as randcplx)
from sympy.abc import x, y, a, b, c, d

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
    assert meijerint_definite(exp(-x), x, 0, x)[0] == \
           (1/x - exp(-x)/x)*exp(I*arg(x))*abs(x)

    # Test -oo to oo
    assert meijerint_definite(exp(-x**2), x, -oo, oo) == (sqrt(pi), True)
    assert meijerint_definite(exp(-abs(x)), x, -oo, oo) == (2, True)
    assert meijerint_definite(exp(-(2*x-3)**2), x, -oo, oo) == (sqrt(pi)/2, True)
    assert meijerint_definite(exp(-abs(2*x-3)), x, -oo, oo) == (1, True)
    assert meijerint_definite(exp(-((x-mu)/sigma)**2/2)/sqrt(2*pi*sigma**2),
                              x, -oo, oo) == (1, True)

def test_lookup_table():
    from random import uniform, randrange
    from sympy import Add
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
            expanded = [hyperexpand(g) for (_, g) in terms]
            assert all (x.is_Piecewise or not x.has(meijerg) for x in expanded)
            expanded = Add(*[f*x for ((f, _), x) in zip(terms, expanded)])
            r = min(abs(formula.subs(subs).n()), abs(expanded.subs(subs).n()))
            if r < 1:
                assert abs(formula.subs(subs).n() - expanded.subs(subs).n()) <= 1e-10
            else:
                assert abs(formula.subs(subs).n() - expanded.subs(subs).n())/r <= 1e-10
