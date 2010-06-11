from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.symbol import Symbol, Wild
from sympy.core.basic import S, C, sympify
from sympy.core.numbers import Rational

from sympy.functions import exp, sin , cos , tan , cot , asin
from sympy.functions import log, sinh, cosh, tanh, coth, asinh

from sympy.functions import sqrt, erf

from sympy.solvers import solve

from sympy.polys import quo, gcd, lcm, \
    monomials, factor, cancel, PolynomialError
from sympy.polys.polyroots import root_factors

from sympy.utilities.iterables import make_list

def components(f, x):
    """Returns a set of all functional components of the given expression
       which includes symbols, function applications and compositions and
       non-integer powers. Fractional powers are collected with with
       minimal, positive exponents.

       >>> from sympy import cos, sin
       >>> from sympy.abc import x, y
       >>> from sympy.integrals.heurisch import components

       >>> components(sin(x)*cos(x)**2, x)
       set([x, cos(x), sin(x)])

    """
    result = set()

    if f.has(x):
        if f.is_Symbol:
            result.add(f)
        elif f.is_Function or f.is_Derivative:
            for g in f.args:
                result |= components(g, x)

            result.add(f)
        elif f.is_Pow:
            result |= components(f.base, x)

            if not f.exp.is_Integer:
                if f.exp.is_Rational:
                    result.add(f.base**Rational(1, f.exp.q))
                else:
                    result |= components(f.exp, x) | set([f])
        else:
            for g in f.args:
                result |= components(g, x)

    return result

# name -> [] of symbols
_symbols_cache = {}

# NB @cacheit is not convenient here
def _symbols(name, n):
    """get vector of symbols local to this module"""
    try:
        lsyms = _symbols_cache[name]
    except KeyError:
        lsyms = []
        _symbols_cache[name] = lsyms

    while len(lsyms) < n:
        lsyms.append( Symbol('%s%i' % (name, len(lsyms)), dummy=True) )

    return lsyms[:n]


def heurisch(f, x, **kwargs):
    """Compute indefinite integral using heuristic Risch algorithm.

       This is a heuristic approach to indefinite integration in finite
       terms using the extended heuristic (parallel) Risch algorithm, based
       on Manuel Bronstein's "Poor Man's Integrator".

       The algorithm supports various classes of functions including
       transcendental elementary or special functions like Airy,
       Bessel, Whittaker and Lambert.

       Note that this algorithm is not a decision procedure. If it isn't
       able to compute the antiderivative for a given function, then this is
       not a proof that such a functions does not exist.  One should use
       recursive Risch algorithm in such case.  It's an open question if
       this algorithm can be made a full decision procedure.

       This is an internal integrator procedure. You should use toplevel
       'integrate' function in most cases,  as this procedure needs some
       preprocessing steps and otherwise may fail.

       Specification
       ============

         heurisch(f, x, rewrite=False, hints=None)

           where
             f : expression
             x : symbol

             rewrite -> force rewrite 'f' in terms of 'tan' and 'tanh'
             hints   -> a list of functions that may appear in anti-derivate

              - hints = None          --> no suggestions at all
              - hints = [ ]           --> try to figure out
              - hints = [f1, ..., fn] --> we know better

       Examples
       ========

       >>> from sympy import tan
       >>> from sympy.integrals.risch import heurisch
       >>> from sympy.abc import x, y

       >>> heurisch(y*tan(x), x)
       y*log(1 + tan(x)**2)/2

       See Manuel Bronstein's "Poor Man's Integrator":

       [1] http://www-sop.inria.fr/cafe/Manuel.Bronstein/pmint/index.html

       For more information on the implemented algorithm refer to:

       [2] K. Geddes, L. Stefanus, On the Risch-Norman Integration
           Method and its Implementation in Maple, Proceedings of
           ISSAC'89, ACM Press, 212-217.

       [3] J. H. Davenport, On the Parallel Risch Algorithm (I),
           Proceedings of EUROCAM'82, LNCS 144, Springer, 144-157.

       [4] J. H. Davenport, On the Parallel Risch Algorithm (III):
           Use of Tangents, SIGSAM Bulletin 16 (1982), 3-6.

       [5] J. H. Davenport, B. M. Trager, On the Parallel Risch
           Algorithm (II), ACM Transactions on Mathematical
           Software 11 (1985), 356-362.

    """
    f = sympify(f)

    if not f.is_Add:
        indep, f = f.as_independent(x)
    else:
        indep = S.One

    if not f.has(x):
        return indep * f * x

    rewritables = {
        (sin, cos, cot)     : tan,
        (sinh, cosh, coth)  : tanh,
    }

    rewrite = kwargs.pop('rewrite', False)

    if rewrite:
        for candidates, rule in rewritables.iteritems():
            f = f.rewrite(candidates, rule)
    else:
        for candidates in rewritables.iterkeys():
            if f.has(*candidates):
                break
        else:
            rewrite = True

    terms = components(f, x)

    hints = kwargs.get('hints', None)

    if hints is not None:
        if not hints:
            a = Wild('a', exclude=[x])
            b = Wild('b', exclude=[x])

            for g in set(terms):
                if g.is_Function:
                    if g.func is exp:
                        M = g.args[0].match(a*x**2)

                        if M is not None:
                            terms.add(erf(sqrt(-M[a])*x))
                elif g.is_Pow:
                    if g.exp.is_Rational and g.exp.q == 2:
                        M = g.base.match(a*x**2 + b)

                        if M is not None and M[b].is_positive:
                            if M[a].is_positive:
                                terms.add(asinh(sqrt(M[a]/M[b])*x))
                            elif M[a].is_negative:
                                terms.add(asin(sqrt(-M[a]/M[b])*x))
        else:
            terms |= set(hints)

    for g in set(terms):
        terms |= components(cancel(g.diff(x)), x)

    V = _symbols('x', len(terms))

    mapping = dict(zip(terms, V))

    rev_mapping = {}

    for k, v in mapping.iteritems():
        rev_mapping[v] = k

    def substitute(expr):
        return expr.subs(mapping)

    diffs = [ substitute(cancel(g.diff(x))) for g in terms ]

    denoms = [ g.as_numer_denom()[1] for g in diffs ]
    try:
        denom = reduce(lambda p, q: lcm(p, q, *V), denoms)
    except PolynomialError:
        # lcm can fail with this. See issue 1418.
        return None

    numers = [ cancel(denom * g) for g in diffs ]

    def derivation(h):
        return Add(*[ d * h.diff(v) for d, v in zip(numers, V) ])

    def deflation(p):
        for y in V:
            if not p.has_any_symbols(y):
                continue

            if derivation(p) is not S.Zero:
                c, q = p.as_poly(y).primitive()
                return deflation(c)*gcd(q, q.diff(y)).as_basic()
        else:
            return p

    def splitter(p):
        for y in V:
            if not p.has_any_symbols(y):
                continue

            if derivation(y) is not S.Zero:
                c, q = p.as_poly(y).primitive()

                q = q.as_basic()

                h = gcd(q, derivation(q), y)
                s = quo(h, gcd(q, q.diff(y), y), y)

                c_split = splitter(c)

                if s.as_poly(y).degree() == 0:
                    return (c_split[0], q * c_split[1])

                q_split = splitter(cancel(q / s))

                return (c_split[0]*q_split[0]*s, c_split[1]*q_split[1])
        else:
            return (S.One, p)

    special = {}

    for term in terms:
        if term.is_Function:
            if term.func is tan:
                special[1 + substitute(term)**2] = False
            elif term.func is tanh:
                special[1 + substitute(term)] = False
                special[1 - substitute(term)] = False
            elif term.func is C.LambertW:
                special[substitute(term)] = True

    F = substitute(f)

    P, Q = F.as_numer_denom()

    u_split = splitter(denom)
    v_split = splitter(Q)

    polys = list(v_split) + [ u_split[0] ] + special.keys()

    s = u_split[0] * Mul(*[ k for k, v in special.iteritems() if v ])
    a, b, c = [ p.as_poly(*V).total_degree() for p in [s, P, Q] ]

    poly_denom = (s * v_split[0] * deflation(v_split[1])).as_basic()

    def exponent(g):
        if g.is_Pow:
            if g.exp.is_Rational and g.exp.q != 1:
                if g.exp.p > 0:
                    return g.exp.p + g.exp.q - 1
                else:
                    return abs(g.exp.p + g.exp.q)
            else:
                return 1
        elif not g.is_Atom:
            return max([ exponent(h) for h in g.args ])
        else:
            return 1

    A, B = exponent(f), a + max(b, c)

    if A > 1 and B > 1:
        monoms = monomials(V, A + B - 1)
    else:
        monoms = monomials(V, A + B)

    poly_coeffs = _symbols('A', len(monoms))

    poly_part = Add(*[ poly_coeffs[i]*monomial
        for i, monomial in enumerate(monoms) ])

    reducibles = set()

    for poly in polys:
        if poly.has(*V):
            try:
                factorization = factor(poly, greedy=True)
            except PolynomialError:
                factorization = poly
            factorization = poly

            if factorization.is_Mul:
                reducibles |= set(factorization.args)
            else:
                reducibles.add(factorization)

    def integrate(field=None):
        irreducibles = set()

        for poly in reducibles:
            for z in poly.atoms(Symbol):
                if z in V:
                    break
            else:
                continue

            irreducibles |= set(root_factors(poly, z, filter=field))

        log_coeffs, log_part = [], []
        B = _symbols('B', len(irreducibles))

        for i, poly in enumerate(irreducibles):
            if poly.has(*V):
                log_coeffs.append(B[i])
                log_part.append(log_coeffs[-1] * log(poly))

        coeffs = poly_coeffs + log_coeffs

        candidate = poly_part/poly_denom + Add(*log_part)

        h = F - derivation(candidate) / denom

        numer = h.as_numer_denom()[0].expand()

        equations = {}

        for term in make_list(numer, Add):
            coeff, dependent = term.as_independent(*V)

            if dependent in equations:
                equations[dependent] += coeff
            else:
                equations[dependent] = coeff

        solution = solve(equations.values(), *coeffs)

        if solution is not None:
            return (solution, candidate, coeffs)
        else:
            return None

    if not (F.atoms(Symbol) - set(V)):
        result = integrate('Q')

        if result is None:
            result = integrate()
    else:
        result = integrate()

    if result is not None:
        (solution, candidate, coeffs) = result

        antideriv = candidate.subs(solution)

        for coeff in coeffs:
            if coeff not in solution:
                antideriv = antideriv.subs(coeff, S.Zero)

        antideriv = antideriv.subs(rev_mapping)
        antideriv = cancel(antideriv).expand()

        if antideriv.is_Add:
            antideriv = antideriv.as_independent(x)[1]

        return indep * antideriv
    else:
        if not rewrite:
            result = heurisch(f, x, rewrite=True, **kwargs)

            if result is not None:
                return indep * result

        return None
