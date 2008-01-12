
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.symbol import Symbol
from sympy.core.function import Function
from sympy.core.basic import Basic, S, Atom
from sympy.core.numbers import Integer, Rational, Zero

from sympy.functions.elementary.trigonometric import sin, cos, cot, tan
from sympy.functions.elementary.hyperbolic import sinh, cosh, tanh, coth

from sympy.solvers import solve
from sympy.simplify.rootof import factors
from sympy.simplify import cancel, simplify, together
from sympy.polynomials import quo, gcd, lcm, factor, PolynomialException

def components(f, x):
    """Returns a set of all functional components of the given expression
       which includes symbols, function applications and compositions and
       non-integer powers. Fractional powers are collected with with
       minimal, positive exponents.

       >>> from sympy import *
       >>> x, y = symbols('xy')

       >>> components(sin(x)*cos(x)**2, x)
       set([sin(x), cos(x), x])

    """
    result = set()

    if f.has(x):
        if isinstance(f, Symbol):
            result.add(f)
        elif isinstance(f, Function):
            for g in f:
                result |= components(g, x)

            result.add(f)
        elif isinstance(f, Pow):
            result |= components(f.base, x)

            if not isinstance(f.exp, Integer):
                if isinstance(f.exp, Rational):
                    result.add(f.base**Rational(1, f.exp.q))
                else:
                    result |= components(f.exp, x) | set([f])
        else:
            for g in f:
                result |= components(g, x)

    return result

def monomials(variables, degree):
    """Generate monomial set of the given total degree or less.

       >>> from sympy import *
       >>> x, y = symbols('xy')

       >>> sorted(monomials([x, y], 2))
       [1, x, y, x**2, y**2, x*y]

       >>> sorted(monomials([x, y], 3))
       [1, x, y, x**2, x**3, y**2, y**3, x*y, x*y**2, y*x**2]

    """
    if not variables:
        return set([S.One])
    else:
        x, tail = variables[0], variables[1:]

        monoms = monomials(tail, degree)

        for i in range(1, degree+1):
            monoms |= set([ x**i * m for m in monomials(tail, degree-i) ])

        return monoms

def symbols(name, n, **kwargs):
    return [ Symbol(name + str(i), **kwargs) for i in range(0, n) ]

def heurisch(f, x, rewrite=False):
    """Compute indefinite integral using heuristic Risch algorithm.

       This is a huristic approach to indefinite integration in finite
       terms using extened heuristic (parallel) Risch algorithm, based
       on Manuel Bronstein's "Poor Man's Integrator".

       The algorithm supports various classes of functions including
       transcendental elementary or special functions like Airy,
       Bessel, Whittaker and Lambert.

       Note that this algorithm is not a decision procedure. If it isn't
       able to compute antiderivative for a given function, then this is
       not a proof that such a functions does not exist. One should use
       recursive Risch algorithm in such case. It's an open question if
       this algorithm can be made a full decision procedure.

       This is an internal integrator procedure. You should use toplevel
       'integrate' function in most cases, as this procedure needs some
       preprocessing steps and otherwise may fail.

       See Manuel Bronstein's "Poor Man's Integrator":

       [1] http://www-sop.inria.fr/cafe/Manuel.Bronstein/pmint/index.html

       For more information on the implemented algorithm refer to:

       [2] K. Geddes, L.Stefanus, On the Risch-Norman Integration
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
    f = Basic.sympify(f)

    indep, f = f.as_independent(x)

    if not f.has(x):
        return indep*f*x

    rewritables = {
        (sin, cos, cot)     : tan,
        (sinh, cosh, coth)  : tanh,
    }

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

    for g in set(terms):
        terms |= components(g.diff(x), x)

    V = symbols('x', len(terms), dummy=True)

    mapping = dict(zip(terms, V))

    rev_mapping = {}

    for k, v in mapping.iteritems():
        rev_mapping[v] = k

    def substitute(expr):
        return expr.subs_dict(mapping)

    diffs = [ substitute(g.diff(x)) for g in terms ]

    denoms = [ g.as_numer_denom()[1] for g in diffs ]
    denom = reduce(lambda p, q: lcm(p, q, V), denoms)

    numers = [ cancel(denom * g, *V) for g in diffs ]

    def derivation(h):
        return Basic.Add(*[ d * h.diff(v) for d, v in zip(numers, V) ])

    def deflation(p):
        for y in p.atoms(Basic.Symbol):
            if not isinstance(derivation(p), Basic.Zero):
                c, q = p.as_polynomial(y).as_primitive()
                return deflation(c) * gcd(q, q.diff(y))
        else:
            return p

    def splitter(p):
        for y in p.atoms(Basic.Symbol):
            if not isinstance(derivation(y), Basic.Zero):
                c, q = p.as_polynomial(y).as_primitive()

                q = q.as_basic()

                h = gcd(q, derivation(q), y)
                s = quo(h, gcd(q, q.diff(y), y), y)

                c_split = splitter(c)

                if s.as_polynomial(y).degree() == 0:
                    return (c_split[0], q * c_split[1])

                q_split = splitter(cancel(q / s, *V))

                return (c_split[0]*q_split[0]*s, c_split[1]*q_split[1])
        else:
            return (S.One, p)

    special = {}

    for term in terms:
        if isinstance(term, Basic.Function):
            if isinstance(term, Basic.tan):
                special[1 + substitute(term)**2] = False
            elif isinstance(term.func, tanh):
                special[1 + substitute(term)] = False
                special[1 - substitute(term)] = False
            #elif isinstance(term.func, Basic.LambertW):
            #    special[substitute(term)] = True

    F = substitute(f)

    P, Q = F.as_numer_denom()

    u_split = splitter(denom)
    v_split = splitter(Q)

    polys = list(v_split) + [ u_split[0] ] + special.keys()

    s = u_split[0] * Basic.Mul(*[ k for k, v in special.iteritems() if v ])
    a, b, c = [ p.as_polynomial(*V).degree() for p in [s, P, Q] ]

    poly_denom = s * v_split[0] * deflation(v_split[1])

    def exponent(g):
        if isinstance(g, Pow):
            if isinstance(g.exp, Rational) and g.exp.q != 1:
                if g.exp.p > 0:
                    return g.exp.p + g.exp.q - 1
                else:
                    return abs(g.exp.p + g.exp.q)
            else:
                return 1
        elif not isinstance(g, Atom):
            return max([ exponent(h) for h in g ])
        else:
            return 1

    A, B = exponent(f), a + max(b, c)

    if A > 1 and B > 1:
        monoms = monomials(V, A + B - 1)
    else:
        monoms = monomials(V, A + B)

    poly_coeffs = symbols('A', len(monoms), dummy=True)

    poly_part = Add(*[ poly_coeffs[i]*monomial
        for i, monomial in enumerate(monoms) ])

    reducibles = set()

    for poly in polys:
        if poly.has(*V):
            try:
                factorization = factor(poly, V)
            except PolynomialException:
                factorization = poly

            if isinstance(factorization, Mul):
                reducibles |= set(factorization[:])
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

            irreducibles |= set(factors(poly,
                z, field, factor=False))

        log_coeffs, log_part = [], []

        for i, poly in enumerate(irreducibles):
            if poly.has(*V):
                log_coeffs.append(Symbol('B%s' % i, dummy=True))
                log_part.append(log_coeffs[-1] * Basic.log(poly))

        coeffs = poly_coeffs + log_coeffs

        candidate = poly_part/poly_denom + Add(*log_part)

        h = together(F - derivation(candidate) / denom)

        numer = h.as_numer_denom()[0].expand()

        if not isinstance(numer, Add):
            numer = [numer]

        equations = {}

        for term in numer:
            coeff, dependent = term.as_independent(*V)

            if dependent in equations:
                equations[dependent] += coeff
            else:
                equations[dependent] = coeff

        solution = solve(equations.values(), coeffs)

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

        antideriv = candidate.subs_dict(solution)

        for coeff in coeffs:
            if coeff not in solution:
                antideriv = antideriv.subs(coeff, S.Zero)

        antideriv = simplify(antideriv.subs_dict(rev_mapping)).expand()

        if isinstance(antideriv, Basic.Add):
            antideriv = Basic.Add(*antideriv.as_coeff_factors()[1])

        return indep * antideriv
    else:
        if not rewrite:
            result = heurisch(f, x, rewrite=True)

            if result is not None:
                return indep * result

        return None
