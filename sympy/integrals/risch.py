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
    monomials, factor, cancel, PolynomialError, Poly
from sympy.polys.polyroots import root_factors

from sympy.utilities.iterables import make_list

def components(f, x):
    """Returns a set of all functional components of the given expression
       which includes symbols, function applications and compositions and
       non-integer powers. Fractional powers are collected with with
       minimal, positive exponents.

       >>> from sympy import cos, sin
       >>> from sympy.abc import x, y
       >>> from sympy.integrals.risch import components

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
    """
    Compute indefinite integral using heuristic Risch algorithm.

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

    [6] M. Bronstein, Symbolic Integration I: Transcendental
       Functions, Second Edition, Springer-Verlag, 2005, pp. 35-70

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
                if g.is_Function and g.func is exp:
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

    V = _symbols('x', len(terms)) # Set of dummy generators to make Poly happy
    One = Poly(1, *V)
    Zero = Poly(0, *V)

    mapping = dict(zip(terms, V))

    rev_mapping = {}

    for k, v in mapping.iteritems():
        rev_mapping[v] = k

    def substitute(expr):
        return expr.subs(mapping)

    diffs = [ substitute(cancel(g.diff(x))) for g in terms ]

    denoms = [ g.as_numer_denom()[1] for g in diffs ]
    denoms = [ Poly(g, *V) for g in denoms ]
    try:
        denom = reduce(lambda p, q: p.lcm(q), denoms)
    except PolynomialError:
        # lcm can fail with this. See issue 1418.
        return None

    numers = [ Poly(denom * g, *V) for g in diffs ]

    def derivation(h):
        return sum([ d * h.diff(v) for d, v in zip(numers, V) ])

    def logderivation(coeff, argument):
        return (coeff*sum([ d * argument.diff(v) for d, v in zip(numers, V)]), argument)

    def deflation(p):
        for y in V:
            if not p.has_any_symbols(y):
                continue

            if derivation(p) != Zero:
                c, q = p.as_poly(y).primitive()
                c, q = Poly(c, *V), Poly(q, *V)
                return deflation(c)*q.gcd(q.diff(y))
        else:
            return p

    def splitter(p):
        for y in V:
            if not p.has_any_symbols(y):
                continue

            if derivation(y) != Zero:
                c, q = p.as_poly(y).primitive()
                c, q = Poly(c, *V), Poly(q, *V)

                h = q.gcd(derivation(q))
                s = h.quo(q.gcd(q.diff(y)))

                c_split = splitter(c)

                if s.degree(y) == 0:
                    return (c_split[0], q * c_split[1])

                q_split = splitter(q.quo(s))

                return (c_split[0]*q_split[0]*s, c_split[1]*q_split[1])
        else:
            return (One, p)

    special = {}

    for term in terms:
        if term.is_Function:
            if term.func is tan:
                special[Poly(1 + substitute(term)**2, *V)] = False
            elif term.func is tanh:
                special[Poly(1 + substitute(term), *V)] = False
                special[Poly(1 - substitute(term), *V)] = False
            elif term.func is C.LambertW:
                special[Poly(substitute(term), *V)] = True

    F = substitute(f)

    P, Q = F.as_numer_denom()
    P, Q = Poly(P, *V), Poly(Q, *V)

    u_split = splitter(denom)
    v_split = splitter(Q)

    polys = list(v_split) + [ u_split[0] ] + special.keys()

    s = u_split[0] * reduce(lambda x, y: x*y,
        [ k for k, v in special.iteritems() if v ], One)
    a, b, c = [ p.total_degree() for p in [s, P, Q] ]

    poly_denom = (s * v_split[0] * deflation(v_split[1]))

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
        monoms = [Poly(t, *V) for t in monomials(V, A + B - 1)]
    else:
        monoms = [Poly(t, *V) for t in monomials(V, A + B)]

    poly_coeffs = _symbols('A', len(monoms))
    all_coeffs = V + poly_coeffs

    poly_part = sum([ Poly(poly_coeffs[i])*monomial.as_poly(*all_coeffs)
        for i, monomial in enumerate(monoms) ])

    reducibles = set()

    for poly in polys:
        if poly.has(*V):

            factorization = poly.factor_list_include() # XXX: greedy=True?

            reducibles |= set([b**e for b, e in factorization])


    def integrate(field=None):
        irreducibles = set()

        for poly in reducibles:
            for z in poly.exclude().gens:
                if z in V:
                    break
            else:
                continue

            irreducibles |= set(root_factors(poly, z, filter=field))

        log_coeffs, log_part, poly_log_part = [], [], []
        B = _symbols('B', len(irreducibles))

        for i, poly in enumerate(irreducibles):
            if poly.has(*V):
                log_coeffs.append(B[i])
                log_part.append(log_coeffs[-1] * log(poly.as_basic()))
                poly_log_part.append(poly.as_poly(*V))

        coeffs = poly_coeffs + log_coeffs

        candidate = poly_part.as_basic()/poly_denom.as_basic() + Add(*log_part)

        ppart_numer_d = poly_denom*derivation(poly_part) - \
            poly_part*derivation(poly_denom)
        ppart_denom_d = poly_denom**2

        logpart_d = [ logderivation(Poly(log_coeffs[i], *coeffs), poly_log_part[i])
            for i in range(len(log_part)) ]
        logpart_denom_d = reduce(lambda a, b: a*b,
            [i for _, i in logpart_d], One)
        logpart_numer_d = [ (logpart_denom_d*i).cancel(logpart_denom_d)
            for i, _ in logpart_d ]
        logpart_numer_d = sum([ Poly(a, *coeffs)*b
            for a, b, _ in logpart_numer_d ])

        candidate_numer_d = ppart_numer_d*logpart_denom_d + \
            logpart_numer_d*ppart_denom_d
        candidate_denom_d = ppart_denom_d*logpart_denom_d

        numer = (P*denom*candidate_denom_d - Q*candidate_numer_d).cancel(denom)
        numer = Poly(numer[0], *coeffs)*numer[1]
        numer = numer.cancel(Q)
        numer = Poly(numer[0], *coeffs)*numer[1]
        numer = numer.cancel(candidate_denom_d)
        numer = Poly(numer[0], *coeffs)*numer[1]


 #       h = F - derivation(candidate) / denom

#        numer = h.as_numer_denom()[0].expand()

        equations = {}

        for term in make_list(numer.as_basic(), Add):
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
        return None
        if not rewrite:
            result = heurisch(f, x, rewrite=True, **kwargs)

            if result is not None:
                return indep * result

        return None
