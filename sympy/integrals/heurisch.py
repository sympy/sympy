from __future__ import print_function, division

from collections import defaultdict
from itertools import permutations

from sympy.core.add import Add
from sympy.core.basic import Basic
from sympy.core.mul import Mul
from sympy.core.symbol import Symbol, Wild, Dummy
from sympy.core.basic import C, sympify
from sympy.core.numbers import Rational, I, pi
from sympy.core.relational import Eq
from sympy.core.singleton import S

from sympy.functions import exp, sin, cos, tan, cot, asin, atan
from sympy.functions import log, sinh, cosh, tanh, coth, asinh, acosh
from sympy.functions import sqrt, erf, erfi, li, Ei
from sympy.functions.elementary.piecewise import Piecewise

from sympy.logic.boolalg import And
from sympy.solvers.solvers import solve, denoms
from sympy.utilities.iterables import uniq

from sympy.polys import quo, gcd, lcm, factor, cancel, PolynomialError
from sympy.polys.monomials import itermonomials
from sympy.polys.polyroots import root_factors

from sympy.polys.rings import PolyRing
from sympy.polys.solvers import solve_lin_sys
from sympy.polys.constructor import construct_domain

from sympy.core.compatibility import reduce, default_sort_key


def components(f, x):
    """
    Returns a set of all functional components of the given expression
    which includes symbols, function applications and compositions and
    non-integer powers. Fractional powers are collected with with
    minimal, positive exponents.

    >>> from sympy import cos, sin
    >>> from sympy.abc import x, y
    >>> from sympy.integrals.heurisch import components

    >>> components(sin(x)*cos(x)**2, x)
    set([x, sin(x), cos(x)])

    See Also
    ========

    heurisch
    """
    result = set()

    if x in f.free_symbols:
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
        lsyms.append( Dummy('%s%i' % (name, len(lsyms))) )

    return lsyms[:n]


def heurisch_wrapper(f, x, rewrite=False, hints=None, mappings=None, retries=3,
                     degree_offset=0, unnecessary_permutations=None):
    """
    A wrapper around the heurisch integration algorithm.

    This method takes the result from heurisch and checks for poles in the
    denominator. For each of these poles, the integral is reevaluated, and
    the final integration result is given in terms of a Piecewise.

    Examples
    ========

    >>> from sympy.core import symbols
    >>> from sympy.functions import cos
    >>> from sympy.integrals.heurisch import heurisch, heurisch_wrapper
    >>> n, x = symbols('n x')
    >>> heurisch(cos(n*x), x)
    sin(n*x)/n
    >>> heurisch_wrapper(cos(n*x), x)
    Piecewise((x, n == 0), (sin(n*x)/n, True))

    See Also
    ========

    heurisch
    """
    f = sympify(f)
    if x not in f.free_symbols:
        return f*x

    res = heurisch(f, x, rewrite, hints, mappings, retries, degree_offset,
                   unnecessary_permutations)
    if not isinstance(res, Basic):
        return res
    # We consider each denominator in the expression, and try to find
    # cases where one or more symbolic denominator might be zero. The
    # conditions for these cases are stored in the list slns.
    slns = []
    for d in denoms(res):
        try:
            slns += solve(d, dict=True, exclude=(x,))
        except NotImplementedError:
            pass
    if not slns:
        return res
    slns = list(uniq(slns))
    # Remove the solutions corresponding to poles in the original expression.
    slns0 = []
    for d in denoms(f):
        try:
            slns0 += solve(d, dict=True, exclude=(x,))
        except NotImplementedError:
            pass
    slns = [s for s in slns if s not in slns0]
    if not slns:
        return res
    if len(slns) > 1:
        eqs = []
        for sub_dict in slns:
            eqs.extend([Eq(key, value) for key, value in sub_dict.items()])
        slns = solve(eqs, dict=True, exclude=(x,)) + slns
    # For each case listed in the list slns, we reevaluate the integral.
    pairs = []
    for sub_dict in slns:
        expr = heurisch(f.subs(sub_dict), x, rewrite, hints, mappings, retries,
                        degree_offset, unnecessary_permutations)
        cond = And(*[Eq(key, value) for key, value in sub_dict.items()])
        pairs.append((expr, cond))
    pairs.append((heurisch(f, x, rewrite, hints, mappings, retries,
                           degree_offset, unnecessary_permutations), True))
    return Piecewise(*pairs)


def heurisch(f, x, rewrite=False, hints=None, mappings=None, retries=3,
             degree_offset=0, unnecessary_permutations=None):
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
    =============

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
    >>> from sympy.integrals.heurisch import heurisch
    >>> from sympy.abc import x, y

    >>> heurisch(y*tan(x), x)
    y*log(tan(x)**2 + 1)/2

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

    See Also
    ========

    sympy.integrals.integrals.Integral.doit
    sympy.integrals.integrals.Integral
    components
    """
    f = sympify(f)
    if x not in f.free_symbols:
        return f*x

    if not f.is_Add:
        indep, f = f.as_independent(x)
    else:
        indep = S.One

    rewritables = {
        (sin, cos, cot): tan,
        (sinh, cosh, coth): tanh,
    }

    if rewrite:
        for candidates, rule in rewritables.items():
            f = f.rewrite(candidates, rule)
    else:
        for candidates in rewritables.keys():
            if f.has(*candidates):
                break
        else:
            rewrite = True

    terms = components(f, x)

    if hints is not None:
        if not hints:
            a = Wild('a', exclude=[x])
            b = Wild('b', exclude=[x])
            c = Wild('c', exclude=[x])

            for g in set(terms):
                if g.is_Function:
                    if g.func is li:
                        M = g.args[0].match(a*x**b)

                        if M is not None:
                            terms.add( x*(li(M[a]*x**M[b]) - (M[a]*x**M[b])**(-1/M[b])*Ei((M[b]+1)*log(M[a]*x**M[b])/M[b])) )
                            #terms.add( x*(li(M[a]*x**M[b]) - (x**M[b])**(-1/M[b])*Ei((M[b]+1)*log(M[a]*x**M[b])/M[b])) )
                            #terms.add( x*(li(M[a]*x**M[b]) - x*Ei((M[b]+1)*log(M[a]*x**M[b])/M[b])) )
                            #terms.add( li(M[a]*x**M[b]) - Ei((M[b]+1)*log(M[a]*x**M[b])/M[b]) )

                    elif g.func is exp:
                        M = g.args[0].match(a*x**2)

                        if M is not None:
                            if M[a].is_positive:
                                terms.add(erfi(sqrt(M[a])*x))
                            else: # M[a].is_negative or unknown
                                terms.add(erf(sqrt(-M[a])*x))

                        M = g.args[0].match(a*x**2 + b*x + c)

                        if M is not None:
                            if M[a].is_positive:
                                terms.add(sqrt(pi/4*(-M[a]))*exp(M[c] - M[b]**2/(4*M[a]))*
                                          erfi(sqrt(M[a])*x + M[b]/(2*sqrt(M[a]))))
                            elif M[a].is_negative:
                                terms.add(sqrt(pi/4*(-M[a]))*exp(M[c] - M[b]**2/(4*M[a]))*
                                          erf(sqrt(-M[a])*x - M[b]/(2*sqrt(-M[a]))))

                        M = g.args[0].match(a*log(x)**2)

                        if M is not None:
                            if M[a].is_positive:
                                terms.add(erfi(sqrt(M[a])*log(x) + 1/(2*sqrt(M[a]))))
                            if M[a].is_negative:
                                terms.add(erf(sqrt(-M[a])*log(x) - 1/(2*sqrt(-M[a]))))

                elif g.is_Pow:
                    if g.exp.is_Rational and g.exp.q == 2:
                        M = g.base.match(a*x**2 + b)

                        if M is not None and M[b].is_positive:
                            if M[a].is_positive:
                                terms.add(asinh(sqrt(M[a]/M[b])*x))
                            elif M[a].is_negative:
                                terms.add(asin(sqrt(-M[a]/M[b])*x))

                        M = g.base.match(a*x**2 - b)

                        if M is not None and M[b].is_positive:
                            if M[a].is_positive:
                                terms.add(acosh(sqrt(M[a]/M[b])*x))
                            elif M[a].is_negative:
                                terms.add((-M[b]/2*sqrt(-M[a])*
                                           atan(sqrt(-M[a])*x/sqrt(M[a]*x**2 - M[b]))))

        else:
            terms |= set(hints)

    for g in set(terms):
        terms |= components(cancel(g.diff(x)), x)

    # TODO: caching is significant factor for why permutations work at all. Change this.
    V = _symbols('x', len(terms))

    mapping = dict(list(zip(terms, V)))

    rev_mapping = {}

    if unnecessary_permutations is None:
        unnecessary_permutations = []
    for k, v in mapping.items():
        rev_mapping[v] = k

    if mappings is None:
        # Pre-sort mapping in order of largest to smallest expressions (last is always x).
        def _sort_key(arg):
            return default_sort_key(arg[0].as_independent(x)[1])
        #optimizing the number of permutations of mappping
        unnecessary_permutations = [(x, mapping[x])]
        del mapping[x]
        mapping = sorted(list(mapping.items()), key=_sort_key, reverse=True)
        mappings = permutations(mapping)

    def _substitute(expr):
        return expr.subs(mapping)

    for mapping in mappings:
        mapping = list(mapping)
        mapping = mapping + unnecessary_permutations
        diffs = [ _substitute(cancel(g.diff(x))) for g in terms ]
        denoms = [ g.as_numer_denom()[1] for g in diffs ]
        if all(h.is_polynomial(*V) for h in denoms) and _substitute(f).is_rational_function(*V):
            denom = reduce(lambda p, q: lcm(p, q, *V), denoms)
            break
    else:
        if not rewrite:
            result = heurisch(f, x, rewrite=True, hints=hints, unnecessary_permutations=unnecessary_permutations)

            if result is not None:
                return indep*result
        return None

    numers = [ cancel(denom*g) for g in diffs ]
    def _derivation(h):
        return Add(*[ d * h.diff(v) for d, v in zip(numers, V) ])

    def _deflation(p):
        for y in V:
            if not p.has(y):
                continue

            if _derivation(p) is not S.Zero:
                c, q = p.as_poly(y).primitive()
                return _deflation(c)*gcd(q, q.diff(y)).as_expr()
        else:
            return p

    def _splitter(p):
        for y in V:
            if not p.has(y):
                continue

            if _derivation(y) is not S.Zero:
                c, q = p.as_poly(y).primitive()

                q = q.as_expr()

                h = gcd(q, _derivation(q), y)
                s = quo(h, gcd(q, q.diff(y), y), y)

                c_split = _splitter(c)

                if s.as_poly(y).degree() == 0:
                    return (c_split[0], q * c_split[1])

                q_split = _splitter(cancel(q / s))

                return (c_split[0]*q_split[0]*s, c_split[1]*q_split[1])
        else:
            return (S.One, p)

    special = {}

    for term in terms:
        if term.is_Function:
            if term.func is tan:
                special[1 + _substitute(term)**2] = False
            elif term.func is tanh:
                special[1 + _substitute(term)] = False
                special[1 - _substitute(term)] = False
            elif term.func is C.LambertW:
                special[_substitute(term)] = True

    F = _substitute(f)

    P, Q = F.as_numer_denom()

    u_split = _splitter(denom)
    v_split = _splitter(Q)

    polys = list(v_split) + [ u_split[0] ] + list(special.keys())

    s = u_split[0] * Mul(*[ k for k, v in special.items() if v ])
    polified = [ p.as_poly(*V) for p in [s, P, Q] ]

    if None in polified:
        return None

    a, b, c = [ p.total_degree() for p in polified ]

    poly_denom = (s * v_split[0] * _deflation(v_split[1])).as_expr()

    def _exponent(g):
        if g.is_Pow:
            if g.exp.is_Rational and g.exp.q != 1:
                if g.exp.p > 0:
                    return g.exp.p + g.exp.q - 1
                else:
                    return abs(g.exp.p + g.exp.q)
            else:
                return 1
        elif not g.is_Atom and g.args:
            return max([ _exponent(h) for h in g.args ])
        else:
            return 1

    A, B = _exponent(f), a + max(b, c)

    if A > 1 and B > 1:
        monoms = itermonomials(V, A + B - 1 + degree_offset)
    else:
        monoms = itermonomials(V, A + B + degree_offset)

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

    def _integrate(field=None):
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

        # TODO: Currently it's better to use symbolic expressions here instead
        # of rational functions, because it's simpler and FracElement doesn't
        # give big speed improvement yet. This is because cancelation is slow
        # due to slow polynomial GCD algorithms. If this gets improved then
        # revise this code.
        candidate = poly_part/poly_denom + Add(*log_part)
        h = F - _derivation(candidate) / denom
        raw_numer = h.as_numer_denom()[0]

        # Rewrite raw_numer as a polynomial in K[coeffs][V] where K is a field
        # that we have to determine. We can't use simply atoms() because log(3),
        # sqrt(y) and similar expressions can appear, leading to non-trivial
        # domains.
        syms = set(coeffs) | set(V)
        non_syms = set([])

        def find_non_syms(expr):
            if expr.is_Integer or expr.is_Rational:
                pass # ignore trivial numbers
            elif expr in syms:
                pass # ignore variables
            elif not expr.has(*syms):
                non_syms.add(expr)
            elif expr.is_Add or expr.is_Mul or expr.is_Pow:
                list(map(find_non_syms, expr.args))
            else:
                # TODO: Non-polynomial expression. This should have been
                # filtered out at an earlier stage.
                raise PolynomialError

        try:
            find_non_syms(raw_numer)
        except PolynomialError:
            return None
        else:
            ground, _ = construct_domain(non_syms, field=True)

        coeff_ring = PolyRing(coeffs, ground)
        ring = PolyRing(V, coeff_ring)

        numer = ring.from_expr(raw_numer)

        solution = solve_lin_sys(numer.coeffs(), coeff_ring)

        if solution is None:
            return None
        else:
            solution = [ (k.as_expr(), v.as_expr()) for k, v in solution.items() ]
            return candidate.subs(solution).subs(list(zip(coeffs, [S.Zero]*len(coeffs))))

    if not (F.atoms(Symbol) - set(V)):
        solution = _integrate('Q')

        if solution is None:
            solution = _integrate()
    else:
        solution = _integrate()

    if solution is not None:
        antideriv = solution.subs(rev_mapping)
        antideriv = cancel(antideriv).expand(force=True)

        if antideriv.is_Add:
            antideriv = antideriv.as_independent(x)[1]

        return indep*antideriv
    else:
        if retries >= 0:
            result = heurisch(f, x, mappings=mappings, rewrite=rewrite, hints=hints, retries=retries - 1, unnecessary_permutations=unnecessary_permutations)

            if result is not None:
                return indep*result

        return None
