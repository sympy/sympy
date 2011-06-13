"""
Integrate functions by rewriting them as Meijer G-functions.

There are three user-visible functions that can be used by other parts of the
sympy library to solve various integration problems:

- meijerint_indefinite
- meijerint_definite
- meijerint_inversion

They can be used to compute, respectively, indefinite integrals, definite
integrals over intervals of the real line, and inverse laplace-type integrals
(from c-I*oo to c+I*oo). See the respective docstrings for details.

The main references for this are:

[L] Luke, Y. L. (1969), The Special Functions and Their Approximations,
    Volume 1

[R] Kelly B. Roach.  Meijer G Function Representations.
    In: Proceedings of the 1997 International Symposium on Symbolic and
    Algebraic Computation, pages 205-211, New York, 1997. ACM.

[P] A. P. Prudnikov, Yu. A. Brychkov and O. I. Marichev (1990).
    Integrals and Series: More Special Functions, Vol. 3,.
    Gordon and Breach Science Publisher
"""
from sympy.core import oo, S, pi
from sympy.core.function import expand, expand_mul
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.symbol import Dummy, Wild
from sympy.simplify import hyperexpand, powdenest
from sympy.logic.boolalg import And, Or
from sympy.functions.special.delta_functions import Heaviside
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.special.hyper import meijerg
from sympy import SYMPY_DEBUG

# keep this at top for easy reference
z = Dummy('z')
def _create_lookup_table(table):
    """ Add formulae for the function -> meijerg lookup table. """
    pass


####################################################################
# First some helper functions.
####################################################################

def _debug(*args):
    if SYMPY_DEBUG:
        for a in args:
            print a,
        print

def _mytype(f, x):
    """ Create a hashable entity describing the type of f. """
    if not f.has(x):
        return ()
    elif f.is_Function:
        return (type(f),)
    else:
        types = [_mytype(a, x) for a in f.args]
        res = []
        for t in types:
            res += list(t)
        res.sort()
        return tuple(res)

# We need to keep track of branches. This argument dummy represent pi.
# We e.g. replace -x by exp(I*argument_dummy)*x. At the end of the day we
# substitute back pi.
argument_dummy = Dummy('pi2', real=True)
def _compute_branch(expr):
    if isinstance(expr, bool):
        return expr
    # TODO this should be avoided where possible
    expr = expr.rewrite('unbranched').subs(argument_dummy, pi)
    if isinstance(expr, bool):
        return expr
    return expr.rewrite('branched', deep=True)

def _add_branch_marker(c):
    """
    Add a 'branch marker' to c.

    >>> from sympy.integrals.meijerint import _add_branch_marker as abm
    >>> from sympy.abc import x, y
    >>> from sympy import I
    >>> abm(x)
    x
    >>> abm(x*y)
    x*y
    >>> abm(-x)
    x*exp(_pi2*I)
    >>> abm(-x*y)
    x*y*exp(_pi2*I)
    >>> abm(I*x)
    x*exp(_pi2*I/2)
    >>> abm(-I*x)
    x*exp(-_pi2*I/2)
    """
    from sympy import arg, exp, I, Mul
    coeff, m = c.as_coeff_mul(*c.free_symbols)
    m = Mul(*m)
    pimult = arg(coeff)/pi
    c = abs(coeff)*exp(I*argument_dummy*pimult)*m
    return c

def _get_coeff_exp(expr, x):
    """
    When expr is known to be of the form c*x**b, with c and/or b possibly 1,
    return c, b.

    >>> from sympy.abc import x, a, b
    >>> from sympy.integrals.meijerint import _get_coeff_exp
    >>> _get_coeff_exp(a*x**b, x)
    (a, b)
    >>> _get_coeff_exp(x, x)
    (1, 1)
    >>> _get_coeff_exp(2*x, x)
    (2, 1)
    >>> _get_coeff_exp(x**3, x)
    (1, 3)
    """
    from sympy import powsimp
    (c, m) = powsimp(expr).as_coeff_mul(x)
    if not m:
        return c, S(0)
    [m] = m
    if m.is_Pow:
        if m.base != x:
            raise ValueError('expr not of form a*x**b')
        return c, m.exp
    elif m == x:
        return c, S(1)
    else:
        raise ValueError('expr not of form a*x**b: %s' % expr)

def _expand_power_base(f):
    """ Expand (3*x)**y to 3**y*x**y. """
    return expand(f, power_base=True, mul=False, power_exp=False,
                     basic=False, multinomial=False, log=False)

def _split_mul(f, x):
    """
    Split expression ``f`` into fac, po, g, where fac is a constant factor,
    po = x**s for some s independent of s, and g is "the rest".

    >>> from sympy.integrals.meijerint import _split_mul
    >>> from sympy import sin
    >>> from sympy.abc import s, x
    >>> _split_mul((3*x)**s*sin(x**2)*x, x)
    (3**s, x*x**s, sin(x**2))
    """
    fac = S(1)
    po  = S(1)
    g   = S(1)
    f = _expand_power_base(f)

    args = Mul.make_args(f)
    for a in args:
        if not a.has(x):
            fac *= a
        elif a.is_Pow and a.base == x:
            po *= a
        elif a == x:
            po *= x
        else:
            g *= a

    return fac, po, g

def _mul_as_two_parts(f):
    """
    Find all the ways to split f into a product of two terms.
    Return None on failure.

    >>> from sympy.integrals.meijerint import _mul_as_two_parts
    >>> from sympy import sin, exp
    >>> from sympy.abc import x
    >>> _mul_as_two_parts(x*sin(x)*exp(x))
    [(x, exp(x)*sin(x)), (exp(x), x*sin(x)), (sin(x), x*exp(x))]
    """
    from sympy.core.compatibility import combinations
    if not f.is_Mul:
        return None
    gs = f.args
    if len(gs) < 2:
        return None

    res = []
    # we now try all ways to split gs into two subsequences
    # TODO this code generates all splittings into two subsequences of equal
    #      length twice!
    for l in range(1, len(gs)//2 + 1):
        for comb in combinations(range(len(gs)), l):
            fac1 = S(1)
            fac2 = S(1)
            for i in range(len(gs)):
                if i in comb:
                    fac1 *= gs[i]
                else:
                    fac2 *= gs[i]
            res += [(fac1, fac2)]

    return res

def _inflate_g(g, n):
    """ Return C, h such that h is a G function of argument z**n and
        g = C*h. """
    # TODO should this be a method of meijerg?
    # See: [L, page 150, equation (5)]
    def inflate(params, n):
        """ (a1, .., ak) -> (a1/n, (a1+1)/n, ..., (ak + n-1)/n) """
        res = []
        for a in params:
            for i in range(n):
                res.append((a + i)/n)
        return res
    v = S(len(g.ap) - len(g.bq))
    C = n**(1 + g.nu + v/2)
    C /= (2*pi)**((n - 1)*g.delta)
    arg = _add_branch_marker(g.argument)**n * n**(n*v)
    return C, meijerg(inflate(g.an, n), inflate(g.aother, n),
                      inflate(g.bm, n), inflate(g.bother, n), arg)

def _flip_g(g):
    """ Turn the G function into one of inverse argument
        (i.e. G(1/x) -> G'(x)) """
    # See [L], section 5.2
    def tr(l): return [1 - a for a in l]
    return meijerg(tr(g.bm), tr(g.bother), tr(g.an), tr(g.aother), 1/g.argument)


def _powdenest(f, force=True):
    """ Like powdenest(f, force=True), but does not turn arg(x) into 0. """
    # XXX should this be in powdenest?
    from sympy import Function, arg
    myarg = Function('arg')
    return powdenest(f.subs(arg, myarg), force=force).subs(myarg, arg)

####################################################################
# Now the "backbone" functions to do actual integration.
####################################################################

def _rewrite_saxena_1(fac, po, g, x):
    """
    Rewrite the integral fac*po*g dx, from zero to infinity, as
    integral fac*G, where G has argument a*x. Note po=x**s.
    Return fac, G.
    """
    _, s = _get_coeff_exp(po, x)
    a, b = _get_coeff_exp(g.argument, x)

    # We substitute t = x**b.
    C = fac/(abs(b)*a**((s+1)/b - 1))
    # Absorb a factor of (at)**((1 + s)/b - 1).
    def tr(l): return [a + (1 + s)/b - 1 for a in l]
    return C, meijerg(tr(g.an), tr(g.aother), tr(g.bm), tr(g.bother), a*x)

def _check_antecedents_1(g, x, helper=False):
    """
    Return a condition under which the mellin transform of g exists.
    Any power of x has already been absorbed into the G function,
    so this is just int_0^\infty g dx.

    See [L, section 5.6.1]. (Note that s=1.)

    If ``helper`` is True, only check if the MT exists at infinity, i.e. if
    int_1^\infty g dx exists.
    """
    from sympy import Eq, Not, ceiling, Ne, re, unbranched_argument as arg
    delta = g.delta
    eta, _ = _get_coeff_exp(g.argument, x)
    m, n, p, q = S([len(g.bm), len(g.an), len(g.ap), len(g.bq)])
    xi    = m + n - p

    if p > q:
        def tr(l): return [1 - x for x in l]
        return _check_antecedents_1(meijerg(tr(g.bm), tr(g.bother),
                                            tr(g.an), tr(g.aother), x/eta),
                                    x)

    tmp = []
    for b in g.bm:
        tmp += [-re(b) < 1]
    for a in g.an:
        tmp += [1 < 1 - re(a)]
    cond_3 = And(*tmp)

    for b in g.bother:
        tmp += [-re(b) < 1]
    for a in g.aother:
        tmp += [1 < 1 - re(a)]
    cond_3_star = And(*tmp)

    cond_4 = (-re(g.nu) + (q + 1 - p)/2 > q - p)

    def debug(*msg):
        _debug(*msg)

    debug('Checking antecedents for 1 function:')
    debug('  delta=%s, eta=%s, m=%s, n=%s, p=%s, q=%s'
           % (delta, eta, m, n, p, q))
    debug('  ap = %s, %s' % (list(g.an), list(g.aother)))
    debug('  bq = %s, %s' % (list(g.bm), list(g.bother)))
    debug('  cond_3=%s, cond_3*=%s, cond_4=%s' % (cond_3, cond_3_star, cond_4))

    conds = []

    # case 1
    case1 = []
    tmp1 = [1 <= n, p < q, 1 <= m]
    tmp2 = [1 <= p, 1 <= m, Eq(q, p + 1), Not(And(Eq(n, 0), Eq(m, p + 1)))]
    tmp3 = [1 <= p, Eq(q, p)]
    for k in range(ceiling(delta/2) + 1):
        tmp3 += [Ne(abs(arg(eta)), (delta - 2*k)*pi)]
    tmp = [delta > 0, abs(arg(eta)) < delta*pi]
    extra = [Ne(eta, 0), cond_3]
    if helper:
        extra = []
    for t in [tmp1, tmp2, tmp3]:
        case1 += [And(*(t + tmp + extra))]
    conds += case1
    debug('  case 1:', case1)

    # case 2
    extra = [cond_3]
    if helper:
        extra = []
    case2 = [And(Eq(n, 0), p + 1 <= m, m <= q,
                  abs(arg(eta)) < delta*pi, *extra)]
    conds += case2
    debug('  case 2:', case2)

    # case 3
    extra = [cond_3, cond_4]
    if helper:
        extra = []
    case3 = [And(p < q, 1 <= m, delta > 0, Eq(abs(arg(eta)), delta*pi), *extra)]
    case3 += [And(p <= q - 2, Eq(delta, 0), Eq(abs(arg(eta)), 0), *extra)]
    conds += case3
    debug('  case 3:', case3)

    # TODO altered cases 4-7

    # extra case from wofram functions site:
    # (reproduced verbatim from prudnikov, section 2.24.2)
    # http://functions.wolfram.com/HypergeometricFunctions/MeijerG/21/02/01/
    case_extra = []
    case_extra += [Eq(p, q), Eq(delta, 0), Eq(arg(eta), 0), Ne(eta, 0)]
    if not helper:
       case_extra += [cond_3]
    s = []
    for a, b in zip(g.ap, g.bq):
        s += [b - a]
    case_extra += [re(Add(*s)) < 0]
    case_extra = And(*case_extra)
    conds += [case_extra]
    debug('  extra case:', [case_extra])

    case_extra_2 = [And(delta > 0, abs(arg(eta)) < delta*pi)]
    if not helper:
        case_extra_2 += [cond_3]
    case_extra_2 = And(*case_extra_2)
    conds += [case_extra_2]
    debug('  second extra case:', [case_extra_2])

    # TODO This leaves only one case from the three listed by prudnikov.
    #      Investigate if these indeed cover everything; if so, remove the rest.

    return Or(*conds)

def _int0oo_1(g, x):
    """
    Evaluate int_0^\infty g dx using G functions,
    assuming the necessary conditions are fulfilled.

    >>> from sympy.abc import a, b, c, d, x, y
    >>> from sympy import meijerg
    >>> from sympy.integrals.meijerint import _int0oo_1
    >>> _int0oo_1(meijerg([a], [b], [c], [d], x*y), x)
    gamma(-a)*gamma(c + 1)/(y*gamma(-d)*gamma(b + 1))
    """
    # See [L, section 5.6.1]. Note that s=1.
    from sympy import gamma, combsimp
    eta, _ = _get_coeff_exp(g.argument, x)
    res = 1/eta
    # XXX TODO we should reduce order first
    for b in g.bm:
        res *= gamma(b + 1)
    for a in g.an:
        res *= gamma(1 - a - 1)
    for b in g.bother:
        res /= gamma(1 - b - 1)
    for a in g.aother:
        res /= gamma(a + 1)
    return combsimp(res)

def _rewrite_saxena(fac, po, g1, g2, x):
    """
    Rewrite the integral fac*po*g1*g2 from 0 to oo in terms of G functions
    with argument c*x.
    Return C, f1, f2 such that integral C f1 f2 from 0 to infinity equals
    integral fac po g1 g2 from 0 to infinity.

    >>> from sympy.integrals.meijerint import _rewrite_saxena
    >>> from sympy.abc import s, t, m
    >>> from sympy import meijerg
    >>> g1 = meijerg([], [], [0], [], s*t)
    >>> g2 = meijerg([], [], [m/2], [-m/2], t**2/4)
    >>> r = _rewrite_saxena(1, t**0, g1, g2, t)
    >>> r[0]
    s/(4*sqrt(pi))
    >>> r[1]
    meijerg(((), ()), ((-1/2, 0), ()), s**2*t/4)
    >>> r[2]
    meijerg(((), ()), ((m/2,), (-m/2,)), t/4)
    """
    from sympy.core.numbers import ilcm
    _, s = _get_coeff_exp(po, x)
    _, b1 = _get_coeff_exp(g1.argument, x)
    _, b2 = _get_coeff_exp(g2.argument, x)
    if b1 < 0:
        b1 = -b1
        g1 = _flip_g(g1)
    if b2 < 0:
        b2 = -b2
        g2 = _flip_g(g2)
    # We may assume b1, b2 are Rationals.
    m1, n1 = b1.p, b1.q
    m2, n2 = b2.p, b2.q
    tau = ilcm(m1*n2, m2*n1)
    r1 = tau//(m1*n2)
    r2 = tau//(m2*n1)

    C1, g1 = _inflate_g(g1, r1)
    C2, g2 = _inflate_g(g2, r2)

    fac *= C1*C2
    a1, b = _get_coeff_exp(g1.argument, x)
    a2, _ = _get_coeff_exp(g2.argument, x)

    # arbitrarily tack on the x**s part to g1
    # TODO should we try both?
    exp = (s + 1)/b - 1
    fac = fac/(abs(b) * a1**exp)
    def tr(l): return [a + exp for a in l]
    g1 = meijerg(tr(g1.an), tr(g1.aother), tr(g1.bm), tr(g1.bother), a1*x)
    g2 = meijerg(g2.an, g2.aother, g2.bm, g2.bother, a2*x)

    return _powdenest(fac, force=True), g1, g2

def _check_antecedents(g1, g2, x):
    """ Return a condition under which the integral theorem applies. """
    from sympy import re, Eq, Not, Ne, cos, I, exp, ceiling, sin, sign
    from sympy import unbranched_argument as arg
    #  Yes, this is madness.
    # XXX TODO this is a testing *nightmare*

    # The following conditions are found in
    # [P], Section 2.24.1
    #
    # They are also reproduced (verbatim!) at
    # http://functions.wolfram.com/HypergeometricFunctions/MeijerG/21/02/03/
    #
    # Note: k=l=r=alpha=1
    sigma, _ = _get_coeff_exp(g1.argument, x)
    omega, _ = _get_coeff_exp(g2.argument, x)
    s, t, u, v = S([len(g1.bm), len(g1.an), len(g1.ap), len(g1.bq)])
    m, n, p, q = S([len(g2.bm), len(g2.an), len(g2.ap), len(g2.bq)])
    bstar = s + t - (u + v)/2
    cstar = m + n - (p + q)/2
    rho = g1.nu + (u - v)/2 + 1
    mu  = g2.nu + (p - q)/2 + 1
    phi = q - p - (v - u)
    eta = 1 - (v - u) - mu - rho
    psi = (pi*(q - m - n) + abs(arg(omega)))/(q - p)
    theta = (pi*(v - s - t) + abs(arg(sigma)))/(v - u)
    lambda_c = (q - p)*abs(omega)**(1/(q - p))*cos(psi) \
             + (v - u)*abs(sigma)**(1/(v - u))*cos(theta)

    def lambda_s0(c1, c2):
        return c1*(q-p)*abs(omega)**(1/(q-p))*sin(psi) \
             + c2*(v-u)*abs(sigma)**(1/(v-u))*sin(theta)
    lambda_s = Piecewise(
        ((lambda_s0(+1, +1)*lambda_s0(-1, -1)),
         And(Eq(arg(sigma), 0), Eq(arg(omega), 0))),
        (lambda_s0(sign(arg(omega)), +1)*lambda_s0(sign(arg(omega)), -1),
         And(Eq(arg(sigma), 0), Ne(arg(omega), 0))),
        (lambda_s0(+1, sign(arg(sigma)))*lambda_s0(-1, sign(arg(sigma))),
         And(Ne(arg(sigma), 0), Eq(arg(omega), 0))),
        (lambda_s0(sign(arg(omega)), sign(arg(sigma))), True))

    _debug('Checking antecedents:')
    _debug('  sigma=%s, s=%s, t=%s, u=%s, v=%s, b*=%s, rho=%s'
           % (sigma, s, t, u, v, bstar, rho))
    _debug('  omega=%s, m=%s, n=%s, p=%s, q=%s, c*=%s, mu=%s,'
           % (omega, m, n, p, q, cstar, mu))
    _debug('  phi=%s, eta=%s, psi=%s, theta=%s' % (phi, eta, psi, theta))

    c1 = True
    for g in [g1, g2]:
        for a in g1.an:
            for b in g1.bm:
                diff = a - b
                if diff > 0 and diff.is_integer:
                    c1 = False

    tmp = []
    for b in g1.bm:
        for d in g2.bm:
            tmp += [re(1 + b + d) > 0]
    c2 = And(*tmp)

    tmp = []
    for a in g1.an:
        for c in g2.an:
            tmp += [re(1 + a + c) < 1 + 1]
    c3 = And(*tmp)

    tmp = []
    for c in g1.an:
        tmp += [(p - q)*re(1 + c - 1) - re(mu) > -S(3)/2]
    c4 = And(*tmp)

    tmp = []
    for d in g1.bm:
        tmp += [(p - q)*re(1 + d) - re(mu) > -S(3)/2]
    c5 = And(*tmp)

    tmp = []
    for c in g2.an:
        tmp += [(u - v)*re(1 + c - 1) - re(rho) > -S(3)/2]
    c6 = And(*tmp)

    tmp = []
    for d in g2.bm:
        tmp += [(u - v)*re(1 + d) - re(rho) > -S(3)/2]
    c7 = And(*tmp)

    c8 = (abs(phi) + 2*re((rho - 1)*(q - p) + (v - u)*(q - p) + (mu - 1)*(v - u)) > 0)
    c9 = (abs(phi) - 2*re((rho - 1)*(q - p) + (v - u)*(q - p) + (mu - 1)*(v - u)) > 0)
    c10 = (abs(arg(sigma)) < bstar*pi)
    c11 = Eq(abs(arg(sigma)), bstar*pi)
    c12 = (abs(arg(omega)) < cstar*pi)
    c13 = Eq(abs(arg(omega)), cstar*pi)

    tmp = []
    z0 = exp(-(bstar + cstar)*pi*I)
    tmp += [Ne(z0*omega/sigma, 1), abs(arg(1 - z0*omega/sigma)) < pi]
    tmp += [Eq(phi, 0), bstar - 1 + cstar <= 0]
    c14 = And(*tmp)

    tmp = []
    z0 = exp(-(bstar + cstar)*pi*I)
    tmp += [Ne(z0*sigma/omega, 1), abs(arg(1 - z0*sigma/omega)) < pi]
    tmp += [Eq(phi, 0), cstar - 1 + bstar <= 0]
    c14_alt = And(*tmp)

    # Since r=k=l=1, in our case there is c14_alt which is the same as calling
    # us with (g1, g2) = (g2, g1). The conditions below enumerate all cases
    # (i.e. we don't have to try arguments reversed by hand), and indeed try
    # all symmetric cases. (i.e. whenever there is a condition involving c14,
    # there is also a dual condition which is exactly what we would get when g1,
    # g2 were interchanged, *but c14 was unaltered*).
    # Hence the following seems correct:
    c14 = Or(c14, c14_alt)

    tmp = [lambda_c > 0,
           And(Eq(lambda_c, 0), Ne(lambda_s, 0), re(eta) > -1),
           And(Eq(lambda_c, 0), Eq(lambda_s, 0), re(eta) > 0)]
    c15 = Or(*tmp)

    for cond, i in [(c1, 1), (c2, 2), (c3, 3), (c4, 4), (c5, 5), (c6, 6),
                    (c7, 7), (c8, 8), (c9, 9), (c10, 10), (c11, 11),
                    (c12, 12), (c13, 13), (c14, 14), (c15, 15)]:
       _debug('  c%s:' % i, cond)


    # We will return Or(*conds)
    conds = []

    def pr(count):
      _debug('  case %s:' % count, conds[-1])
    conds += [And(m*n*s*t != 0, bstar > 0, cstar > 0, c1, c2, c3, c10, c12)] #1
    pr(1)
    conds += [And(Eq(u, v), Eq(bstar, 0), cstar > 0, sigma > 0, re(rho) < 1,
                  c1, c2, c3, c12)] #2
    pr(2)
    conds += [And(Eq(p, q), Eq(cstar, 0), bstar > 0, omega > 0, re(mu) < 1,
                  c1, c2, c3, c10)] #3
    pr(3)
    conds += [And(Eq(p, q), Eq(u, v), Eq(bstar, 0), Eq(cstar, 0),
                  sigma > 0, omega > 0, re(mu) < 1, re(rho) < 1,
                  Ne(sigma, omega), c1, c2, c3)] #4
    pr(4)
    conds += [And(Eq(p, q), Eq(u, v), Eq(bstar, 0), Eq(cstar, 0),
                  sigma > 0, omega > 0, re(mu + rho) < 1,
                  Ne(omega, sigma), c1, c2, c3)] #5
    pr(5)
    conds += [And(p > q, s > 0, bstar > 0, cstar >= 0,
                  c1, c2, c3, c5, c10, c13)] #6
    pr(6)
    conds += [And(p < q, t > 0, bstar > 0, cstar >= 0,
                  c1, c2, c3, c4, c10, c13)] #7
    pr(7)
    conds += [And(u > v, m > 0, cstar > 0, bstar >= 0,
                  c1, c2, c3, c7, c11, c12)] #8
    pr(8)
    conds += [And(u < v, n > 0, cstar > 0, bstar >= 0,
                  c1, c2, c3, c6, c11, c12)] #9
    pr(9)
    conds += [And(p > q, Eq(u, v), Eq(bstar, 0), cstar >= 0, sigma > 0,
                  re(rho) < 1, c1, c2, c3, c5, c13)] #10
    pr(10)
    conds += [And(p < q, Eq(u, v), Eq(bstar, 0), cstar >= 0, sigma > 0,
                  re(rho) < 1, c1, c2, c3, c4, c13)] #11
    pr(11)
    conds += [And(Eq(p, q), u > v, bstar >= 0, Eq(cstar, 0), omega > 0,
                  re(mu) < 1, c1, c2, c3, c7, c11)] #12
    pr(12)
    conds += [And(Eq(p, q), u < v, bstar >= 0, Eq(cstar, 0), omega > 0,
                  re(mu) < 1, c1, c2, c3, c6, c11)] #13
    pr(13)
    conds += [And(p < q, u > v, bstar >= 0, cstar >= 0,
                  c1, c2,c3, c4, c7, c11, c13)] #14
    pr(14)
    conds += [And(p > q, u < v, bstar >= 0, cstar >= 0,
                  c1, c2, c3, c5, c6, c11, c13)] #15
    pr(15)
    conds += [And(p > q, u > v, bstar >= 0, cstar >= 0,
                  c1, c2, c3, c5, c7, c8, c11, c13, c14)] #16
    pr(16)
    conds += [And(p < q, u < v, bstar >= 0, cstar >= 0,
                  c1, c2, c3, c4, c6, c9, c11, c13, c14)] #17
    pr(17)
    conds += [And(Eq(t, 0), s > 0, bstar > 0, phi > 0, c1, c2, c10)] #18
    pr(18)
    conds += [And(Eq(s, 0), t > 0, bstar > 0, phi < 0, c1, c3, c10)] #19
    pr(19)
    conds += [And(Eq(n, 0), m > 0, cstar > 0, phi < 0, c1, c2, c12)] #20
    pr(20)
    conds += [And(Eq(m, 0), n > 0, cstar > 0, phi > 0, c1, c3, c12)] #21
    pr(21)
    conds += [And(Eq(s*t, 0), bstar > 0, cstar > 0,
                  c1, c2, c3, c10, c12)] #22
    pr(22)
    conds += [And(Eq(m*n, 0), bstar > 0, cstar > 0,
                  c1, c2, c3, c10, c12)] #23
    pr(23)

    # Let's short-circuit if this worked ...
    # the rest is corner-cases and terrible to read.
    r = Or(*conds)
    if _compute_branch(r) is not False:
        return r

    conds += [And(m + n > p, Eq(t, 0), Eq(phi, 0), s > 0, bstar > 0, cstar < 0,
                  abs(arg(omega)) < (m + n - p + 1)*pi,
                  c1, c2, c10, c14, c15)] #24
    pr(24)
    conds += [And(m + n > q, Eq(s, 0), Eq(phi, 0), t > 0, bstar > 0, cstar < 0,
                  abs(arg(omega)) < (m + n - q + 1)*pi,
                  c1, c3, c10, c14, c15)] #25
    pr(25)
    conds += [And(Eq(p, q - 1), Eq(t, 0), Eq(phi, 0), s > 0, bstar > 0,
                  cstar >= 0, cstar*pi < abs(arg(omega)),
                  c1, c2, c10, c14, c15)] #26
    pr(26)
    conds += [And(Eq(p, q + 1), Eq(s, 0), Eq(phi, 0), t > 0, bstar > 0,
                  cstar >= 0, cstar*pi < abs(arg(omega)),
                  c1, c3, c10, c14, c15)] #27
    pr(27)
    conds += [And(p < q - 1, Eq(t, 0), Eq(phi, 0), s > 0, bstar > 0,
                  cstar >= 0, cstar*pi < abs(arg(omega)),
                  abs(arg(omega)) < (m + n - p + 1)*pi,
                  c1, c2, c10, c14, c15)] #28
    pr(28)
    conds += [And(p > q + 1, Eq(s, 0), Eq(phi, 0), t > 0, bstar > 0, cstar >= 0,
                  cstar*pi < abs(arg(omega)),
                  abs(arg(omega)) < (m + n - q + 1)*pi,
                  c1, c3, c10, c14, c15)] #29
    pr(29)
    conds += [And(Eq(n, 0), Eq(phi, 0), s + t > 0, m > 0, cstar > 0, bstar < 0,
                  abs(arg(sigma)) < (s + t - u + 1)*pi,
                  c1, c2, c12, c14, c15)] #30
    pr(30)
    conds += [And(Eq(m, 0), Eq(phi, 0), s + t > v, n > 0, cstar > 0, bstar < 0,
                  abs(arg(sigma)) < (s + t - v + 1)*pi,
                  c1, c3, c12, c14, c15)] #31
    pr(31)
    conds += [And(Eq(n, 0), Eq(phi, 0), Eq(u, v - 1), m > 0, cstar > 0,
                  bstar >= 0, bstar*pi < abs(arg(sigma)),
                  abs(arg(sigma)) < (bstar + 1)*pi,
                  c1, c2, c12, c14, c15)] #32
    pr(32)
    conds += [And(Eq(m, 0), Eq(phi, 0), Eq(u, v + 1), n > 0, cstar > 0,
                  bstar >= 0, bstar*pi < abs(arg(sigma)),
                  abs(arg(sigma)) < (bstar + 1)*pi,
                  c1, c3, c12, c14, c15)] #33
    pr(33)
    conds += [And(Eq(n, 0), Eq(phi, 0), u < v - 1, m > 0, cstar > 0, bstar >= 0,
                  bstar*pi < abs(arg(sigma)),
                  abs(arg(sigma)) < (s + t - u + 1)*pi,
                  c1, c2, c12, c14, c15)] #34
    pr(34)
    conds += [And(Eq(m, 0), Eq(phi, 0), u > v + 1, n > 0, cstar > 0, bstar >= 0,
                  bstar*pi < abs(arg(sigma)),
                  abs(arg(sigma)) < (s + t - v + 1)*pi,
                  c1, c3, c12, c14, c15)] #35
    pr(35)

    # The following case is from [Luke1969]. As far as I can tell, it is *not*
    # covered by prudnikov's.
    # Let G1 and G2 be the two G-functions. Suppose the integral exists from
    # 0 to a > 0 (this is easy the easy part), that G1 is exponential decay at
    # infinity, and that the mellin transform of G2 exists.
    # Then the integral exists.
    mt1_exists = _check_antecedents_1(g1, x, helper=True)
    mt2_exists = _check_antecedents_1(g2, x, helper=True)
    conds += [And(mt2_exists, Eq(t, 0), u < s, bstar > 0, c10, c1, c2, c3)]
    pr('E1')
    conds += [And(mt2_exists, Eq(s, 0), v < t, bstar > 0, c10, c1, c2, c3)]
    pr('E2')
    conds += [And(mt1_exists, Eq(n, 0), p < m, cstar > 0, c12, c1, c2, c3)]
    pr('E3')
    conds += [And(mt1_exists, Eq(m, 0), q < n, cstar > 0, c12, c1, c2, c3)]
    pr('E4')

    return Or(*conds)

    # NOTE An alternative, but as far as I can tell weaker, set of conditions
    #      can be found in [L, section 5.6.2].

def _int0oo(g1, g2, x):
    """
    Express integral from zero to infinity g1*g2 using a G function,
    assuming the necessary conditions are fulfilled.

    >>> from sympy.integrals.meijerint import _int0oo
    >>> from sympy.abc import s, t, m
    >>> from sympy import meijerg, S
    >>> g1 = meijerg([], [], [-S(1)/2, 0], [], s**2*t/4)
    >>> g2 = meijerg([], [], [m/2], [-m/2], t/4)
    >>> _int0oo(g1, g2, t)
    4*meijerg(((1/2, 0), ()), ((m/2,), (-m/2,)), s**(-2))/s**2
    """
    # See: [L, section 5.6.2, equation (1)]
    eta, _   = _get_coeff_exp(g1.argument, x)
    omega, _ = _get_coeff_exp(g2.argument, x)
    def neg(l): return [-x for x in l]
    a1 = neg(g1.bm) + list(g2.an)
    a2 = list(g2.aother) + neg(g1.bother)
    b1 = neg(g1.an) + list(g2.bm)
    b2 = list(g2.bother) + neg(g1.aother)
    return meijerg(a1, a2, b1, b2, omega/eta)/eta


####################################################################
# Finally, the real meat.
####################################################################

_lookup_table = None
def _rewrite_single(f, x):
    """
    Try to rewrite f as a sum of single G functions of the form
    C*x**s*G(a*x**b), where b is a rational number and C is independent of x.
    We guarantee that result.argument.as_coeff_mul(x) returns (a, (x**b,))
    or (a, ()).
    Returns a list of tuples (C, s, G) and a condition cond.
    Returns None on failure.
    """
    global _lookup_table
    if not _lookup_table:
        _lookup_table = {}
        _create_lookup_table(_lookup_table)

    if isinstance(f, meijerg):
        from sympy import factor
        coeff, m = factor(f.argument, x).as_coeff_mul(x)
        if len(m) > 1:
            return None
        m = m[0]
        if m.is_Pow:
            if m.base != x or not m.exp.is_Rational:
                return None
        elif m != x:
            return None
        return [(1, 0, meijerg(f.an, f.aother, f.bm, f.bother, coeff*m))], True

    from sympy.core.relational import Relational
    def tr(expr):
        if expr.rel_op == '<':
            return expr.lhs < expr.rhs
        else:
            raise NotImplementedError

    f = f.subs(x, z)
    t = _mytype(f, z)
    if t in _lookup_table:
        l = _lookup_table[t]
        for formula, terms, cond in l:
            subs = f.match(formula)
            if subs:
                if not isinstance(cond, bool):
                    cond = cond.subs(subs)
                if not isinstance(cond, bool):
                    # XXX is there a better way?
                    cond = cond.replace(lambda x: x.is_Relational, tr)
                if cond is False:
                    continue
                if not isinstance(terms, list):
                    terms = terms(subs)
                return [(fac.subs(subs), 0, g.subs(subs).subs(z, x))
                        for (fac, g) in terms], cond

    # TODO recursive mellin transform

def _rewrite1(f, x):
    """
    Try to rewrite f using a (sum of) single G functions with argument a*x**b.
    Return fac, po, g such that f = fac*po*g, fac is independent of x
    and po = x**s.
    Here g is a result from _rewrite_single.
    Return None on failure.
    """
    fac, po, g = _split_mul(f, x)
    g = _rewrite_single(g, x)
    if g:
        return fac, po, g[0], g[1]

def _rewrite2(f, x):
    """
    Try to rewrite f as a product of two G functions of arguments a*x**b.
    Return fac, po, g1, g2 such that f = fac*po*g1*g2, where fac is
    independent of x and po is x**s.
    Here g1 and g2 are results of _rewrite_single.
    Returns None on failure.
    """
    fac, po, g = _split_mul(f, x)
    l = _mul_as_two_parts(g)
    if not l:
        return None
    for fac1, fac2 in l:
        g1 = _rewrite_single(fac1, x)
        g2 = _rewrite_single(fac2, x)
        if g1 and g2:
            cond = And(g1[1], g2[1])
            if cond is not False:
                return fac, po, g1[0], g2[0], cond


def meijerint_indefinite(f, x):
    """
    Compute an indefinite integral of ``f`` by rewriting it as a G function.
    """
    from sympy import Integral
    _debug('Trying to compute the indefinite integral of', f, 'wrt', x)
    gs = _rewrite1(f, x)
    if gs is None:
        # Note: the code that calls us will do expand() and try again
        return None

    fac, po, gl, cond = gs
    _debug(' could rewrite:', gs)
    res = S(0)
    for C, s, g in gl:
        a, b = _get_coeff_exp(g.argument, x)
        _, c = _get_coeff_exp(po, x)
        c += s

        # we do a substitution t=a*x**b, get integrand fac*t**rho*g
        fac_ = fac * C / (b*a**((1 + c)/b))
        rho = (c + 1)/b - 1

        # we now use t**rho*G(params, t) = G(params + rho, t)
        # [L, page 150, equation (4)]
        # and integral G(params, t) dt = G(1, params+1, 0, t)
        #   (or a similar expression with 1 and 0 exchanged ... pick the one
        #    which yields a well-defined function)
        # [R, section 5]
        t = Dummy('t')
        def tr(p): return [a + rho + 1 for a in p]
        if any(b.is_integer and b <= 0 for b in tr(g.bm)):
            r = -meijerg(tr(g.an), tr(g.aother) + [1], tr(g.bm) + [0], tr(g.bother), t)
        else:
            r = meijerg(tr(g.an) + [1], tr(g.aother), tr(g.bm), tr(g.bother) + [0], t)
        r = hyperexpand(r)

        # now substitute back
        # Note: we really do want to powers of x to combine.
        res += fac_*_powdenest(r.subs(t, a*x**b), force=True)

    # This multiplies out superfluous powers of x we created, and chops off
    # constants (e.g. x*(exp(x)/x - 1/x) -> exp(x))
    res = _compute_branch(Add(*expand_mul(res).as_coeff_add(x)[1]))
    return Piecewise((res, _compute_branch(cond)), (Integral(f, x), True))

def meijerint_definite(f, x, a, b):
    """
    Integrate ``f`` over the interval [``a``, ``b``], by rewriting it as a product
    of two G functions, or as a single G function.

    Return res, cond, where cond are convergence conditions.

    This function is implemented as a succession of functions
    meijerint_definite, _meijerint_definite_2, _meijerint_definite_3,
    _meijerint_definite_4. Each function in the list calls the next one
    (presumably) several times. This means that calling meijerint_definite
    can be very costly.
    """
    # This consists of three steps:
    # 1) Change the integration limits to 0, oo
    # 2) Rewrite in terms of G functions
    # 3) Evaluate the integral
    #
    # There are usually several ways of doing this, and we want to try all.
    # This function does (1), calls _meijerint_definite_2 for step (2).
    from sympy import Integral, arg, exp, I, Wild, And, DiracDelta
    _debug('Integrating', f, 'wrt %s from %s to %s.' % (x, a, b))

    if f.has(DiracDelta):
        _debug('Integrand has DiracDelta terms - giving up.')
        return None

    f_, x_, a_, b_ = f, x, a, b

    # Let's use a dummy in case any of the boundaries has x.
    d = Dummy('x')
    f = f.subs(x, d)
    x = d

    if a == -oo and b != oo:
        return meijerint_definite(f.subs(x, -x), x, -b, -a)

    if a == -oo:
        # Integrating -oo to oo. We need to find a place to split the integral.
        _debug('  Integrating -oo to +oo.')
        p, q = map(lambda n: Wild(n, exclude=[x]), 'pq')
        def compute_innermost(expr, res):
            m = expr.match(p*x+q)
            if m and m[p] != 0:
                res.add(-m[q]/m[p])
                return
            if expr.is_Atom:
                return
            for arg in expr.args:
                compute_innermost(arg, res)
        innermost = set()
        compute_innermost(f, innermost)
        _debug('  Sensible splitting points:', innermost)
        for c in list(innermost) + [S(0)]:
            _debug('  Trying to split at', c)
            if not c.is_real:
                _debug('  Non-real splitting point.')
                continue
            res1 = _meijerint_definite_2(f.subs(x, x + c), x)
            if res1 is None:
                _debug('  But could not compute first integral.')
                continue
            res2 = _meijerint_definite_2(f.subs(x, c-x), x)
            if res2 is None:
                _debug('  But could not compute second integral.')
                continue
            res1, cond1 = res1
            res2, cond2 = res2
            cond = _compute_branch(And(cond1, cond2))
            if cond is False:
                _debug('  But combined condition is always false.')
                continue
            res = _compute_branch(hyperexpand(res1 + res2))
            return res, cond
        return

    f = f.subs(x, x + a)
    b = b - a
    a = 0
    if b != oo:
        phi = exp(I*arg(b))
        b = abs(b)
        f = f.subs(x, phi*x)
        f *= Heaviside(b - x)*phi
        b = oo

    _debug('Changed limits to', a, b)
    _debug('Changed function to', f)
    res = _meijerint_definite_2(f, x)
    if res is not None:
        res, cond = res
        return _compute_branch(hyperexpand(res)), cond

def _guess_expansion(f, x):
    """ Try to guess sensible rewritings for integrand f(x). """
    from sympy import expand_trig
    from sympy.functions.elementary.trigonometric import TrigonometricFunction
    from sympy.functions.elementary.hyperbolic import HyperbolicFunction
    res = [(f, 'originial integrand')]

    expanded = expand_mul(res[-1][0])
    if expanded != res[-1][0]:
        res += [(expanded, 'expand_mul')]

    expanded = expand(res[-1][0])
    if expanded != res[-1][0]:
        res += [(expanded, 'expand')]

    if res[-1][0].has(TrigonometricFunction, HyperbolicFunction):
        expanded = expand_mul(expand_trig(res[-1][0]))
        if expanded != res[-1][0]:
            res += [(expanded, 'expand_trig, expand_mul')]

    return res

def _meijerint_definite_2(f, x):
    """
    Try to integrate f dx from zero to infinty.

    The body of this function computes various 'simplifications'
    f1, f2, ... of f (e.g. by calling expand_mul(), trigexpand()
    - see _guess_expansion) and calls _meijerint_definite_3 with each of
    these in succession.
    If _meijerint_definite_3 succeedes with any of the simplified functions,
    returns this result.
    """
    # This function does preparation for (2), calls
    # _meijerint_definite_3 for (2) and (3) combined.

    # use a positive dummy - we integrate from 0 to oo
    dummy = Dummy('x', positive=True)
    f = f.subs(x, dummy)
    x = dummy

    if f == 0:
        return 0, True

    for g, explanation in _guess_expansion(f, x):
        _debug('Trying', explanation)
        res = _meijerint_definite_3(g, x)
        if res is not None and res[1] is not False:
            return res

def _meijerint_definite_3(f, x):
    """
    Try to integrate f dx from zero to infinity.

    This function calls _meijerint_definite_4 to try to compute the
    integral. If this fails, it tries using linearity.
    """
    res = _meijerint_definite_4(f, x)
    if res is not None and res[1] is not False:
        return res
    if f.is_Add:
        _debug('Expanding and evaluating all terms.')
        ress = [_meijerint_definite_4(g, x) for g in f.args]
        if all(r is not None for r in ress):
            conds = []
            res = S(0)
            for r, c in ress:
                res += r
                conds += [c]
            c = And(*conds)
            if c is not False:
                return res, c

def _meijerint_definite_4(f, x, only_double=False):
    """
    Try to integrate f dx from zero to infinity.

    This function tries to apply the integration theorems found in literature,
    i.e. it tries to rewrite f as either one or a product of two G-functions.

    The parameter ``only_double`` is used internally in the recursive algorithm
    to disable trying to rewrite f as a single G-function.
    """
    # This function does (2) and (3)
    _debug('Integrating', f)
    # Try single G function.
    if not only_double:
        gs = _rewrite1(f, x)
        if gs is not None:
            fac, po, g, cond = gs
            _debug('Could rewrite as single G function:', fac, po, g)
            res = S(0)
            for C, s, f in g:
                if C == 0:
                    continue
                C, f = _rewrite_saxena_1(fac*C, po*x**s, f, x)
                res += C*_int0oo_1(f, x)
                cond = And(cond, _check_antecedents_1(f, x))
            cond = _compute_branch(cond)
            if cond is False:
                _debug('But cond is always False.')
            else:
                return res, cond

    # Try two G functions.
    gs = _rewrite2(f, x)
    if gs is not None:
        fac, po, g1, g2, cond = gs
        _debug('Could rewrite as two G functions:', fac, po, g1, g2)
        res = S(0)
        for C1, s1, f1 in g1:
            for C2, s2, f2 in g2:
                C, f1_, f2_ = _rewrite_saxena(fac*C1*C2, po*x**(s1 + s2), f1, f2, x)
                _debug('Saxena subst for yielded:', C, f1_, f2_)
                cond = And(cond, _check_antecedents(f1_, f2_, x))
                res += C*_int0oo(f1_, f2_, x)
        _debug('Result before branch substitutions is:', res)
        cond = _compute_branch(cond)
        if cond is False:
            _debug('But cond is always False.')
        else:
            return res, cond
