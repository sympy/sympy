
from sympy.core import Basic, S, Symbol


from sympy.solvers import solve
from sympy.polynomials import quo, gcd, lcm, roots, factor_
from sympy.simplify import normal, simplify, trigsimp, together

def components(expr):
    """Returns a set of all functional components of the given expression
       with symbols, functions applications and compositions. All integer
       powers are being skipped however fractional and functional powers
       are collected as well.

       >>> from sympy import *
       >>> x, y = symbols('xy')

       >>> components(x*y)
       set([y, x])

       >>> components(sin(x))
       set([sin(x), x])

       >>> components(sin(x)*cos(x)**2)
       set([sin(x), cos(x), x])

       >>> components(sin(x)*sqrt(log(x)))
       set([log(x), sin(x), log(x)**(1/2), x])

       >>> components(x*sin(exp(x)*y))
       set([y, sin(y*exp(x)), x, exp(x)])

    """
    result = set()

    from sympy.core.function import Function
    if isinstance(expr, Basic.Symbol):
        result.add(expr)
    elif isinstance(expr, Function):
        for obj in expr:
            result |= components(obj)

        result.add(expr)
    elif isinstance(expr, Basic.Pow):
        result |= components(expr.base)

        if not isinstance(expr.exp, Basic.Integer):
            result |= components(expr.exp) | set([expr])
    else:
        for obj in expr:
            result |= components(obj)

    return result

def monomials(variables, degree):
    """Generate monomials set of the given total degree or less.

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

def factorization(poly, linear=False):
    """Returns a list with polynomial factors over rationals or, if
       'linear' flag is set over Q(i). This simple handler should
       be merged with original factor() method.

       >>> from sympy import *
       >>> x, y = symbols('xy')

       >>> factorization(x**2 - y**2)
       set([1, x - y, x + y])

       >>> factorization(x**2 + 1)
       set([1, 1 + x**2])

       >>> factorization(x**2 + 1, linear=True)
       set([1, x - I, I + x])

    """
    factored = set([ q.as_basic() for q in factor_.factor(poly) ])

    if not linear:
        return factored
    else:
        factors = []

        for factor in factored:
            symbols = factor.atoms(Symbol)

            if len(symbols) == 1:
                x = symbols.pop()
            else:
                factors += [ factor ]
                continue

            linearities = [ x - r for r in roots(factor, x) ]

            if not linearities:
                factors += [ factor ]
            else:
                unfactorable = quo(factor, Basic.Mul(*linearities))

                if not isinstance(unfactorable, Basic.Number):
                    factors += [ unfactorable ]

                factors += linearities

        return set(factors)

def risch_norman(f, x, rewrite=False):
    """Computes indefinite integral using extended Risch-Norman algorithm,
       also known as parallel Risch. This is a simplified version of full
       recursive Risch algorithm. It is designed for integrating various
       classes of functions including transcendental elementary or special
       functions like Airy, Bessel, Whittaker and Lambert.

       The main difference between this algorithm and the recursive one
       is that rather than computing a tower of differential extensions
       in a recursive way, it handles all cases in one shot. That's why
       it is called parallel Risch algorithm. This makes it much faster
       than the original approach.

       Another benefit is that it doesn't require to rewrite expressions
       in terms of complex exponentials. Rather it uses tangents and so
       antiderivatives are being found in a more familliar form.

       Risch-Norman algorithm can also handle special functions very
       easily without any additional effort. Just differentiation
       method must be known for a given function.

       Note that this algorithm is not a decision procedure. If it
       computes an antiderivative for a given integral then it's a
       proof that such function exists. However when it fails then
       there still may exist an antiderivative and a fallback to
       recurrsive Risch algorithm would be necessary.

       The question if this algorithm can be made a full featured
       decision procedure still remains open.

       For more information on the implemented algorithm refer to:

       [1] K. Geddes, L.Stefanus, On the Risch-Norman Integration
           Method and its Implementation in Maple, Proceedings of
           ISSAC'89, ACM Press, 212-217.

       [2] J. H. Davenport, On the Parallel Risch Algorithm (I),
           Proceedings of EUROCAM'82, LNCS 144, Springer, 144-157.

       [3] J. H. Davenport, On the Parallel Risch Algorithm (III):
           Use of Tangents, SIGSAM Bulletin 16 (1982), 3-6.

       [4] J. H. Davenport, B. M. Trager, On the Parallel Risch
           Algorithm (II), ACM Transactions on Mathematical
           Software 11 (1985), 356-362.

    """
    f = Basic.sympify(f)

    from sympy import sin, cos, cot, tan, \
                      sinh, cosh, tanh, coth

    if not f.has(x):
        return f * x

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

    terms = components(f)

    for g in set(terms):
        h = g.diff(x)

        if not isinstance(h, Basic.Zero):
            terms |= components(h)

    terms = [ g for g in terms if g.has(x) ]
    #XXX: This fixes bugs: it seems integration variable must be to the right
    #i.e. in terms [exp(x), x], but not [x, exp(x)], this should be checked
    #and made more robust.
    # simply move integration variable to the end of terms
    # FIXME: this is probably incorrect! (--kirr)
    terms.remove(x)
    terms.append(x)

    V, in_terms, out_terms = [], [], {}

    for i, term in enumerate(terms):
        V += [ Symbol('x%s' % i) ]

        N = term.count_ops(symbolic=False)
        in_terms += [ (N, term, V[-1]) ]

        out_terms[V[-1]] = term

    in_terms.sort(lambda u, v: int(v[0] - u[0]))

    def substitute(expr):
        for _, g, symbol in in_terms:
            expr = expr.subs(g, symbol)

        return expr

    diffs = [ substitute(g.diff(x)) for g in terms ]

    denoms = [ g.as_numer_denom()[1] for g in diffs ]
    denom = reduce(lambda p, q: lcm(p, q, V), denoms)

    numers = [ normal(denom * g, *V) for g in diffs ]

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

                q_split = splitter(normal(q / s, *V))

                return (c_split[0]*q_split[0]*s, c_split[1]*q_split[1])
        else:
            return (S.One, p)

    special = []

    for term in terms:
        if isinstance(term, Basic.Function):
            if isinstance(term, Basic.tan):
                special += [ (1 + substitute(term)**2, False) ]
            elif isinstance(term.func, tanh):
                special += [ (1 + substitute(term), False),
                             (1 - substitute(term), False) ]
            #elif isinstance(term.func, Basic.LambertW):
            #    special += [ (substitute(term), True) ]

    ff = substitute(f)

    P, Q = ff.as_numer_denom()

    u_split = splitter(denom)
    v_split = splitter(Q)

    s = u_split[0] * Basic.Mul(*[ g for g, a in special if a ])
    a, b, c = [ p.as_polynomial(*V).degree() for p in [s, P, Q] ]

    candidate_denom = s * v_split[0] * deflation(v_split[1])
    monoms = monomials(V, 1 + a + max(b, c))

    linear = False

    while True:
        coeffs, candidate, factors = [], S.Zero, set()

        for i, monomial in enumerate(monoms):
            coeffs += [ Symbol('A%s' % i) ]
            candidate += coeffs[-1] * monomial

        candidate /= candidate_denom

        polys = [ v_split[0], v_split[1], u_split[0]] + [ s[0] for s in special ]

        for irreducibles in [ factorization(p, linear) for p in polys ]:
            factors |= irreducibles

        for i, irreducible in enumerate(factors):
            if not isinstance(irreducible, Basic.Number):
                coeffs += [ Symbol('B%s' % i) ]
                candidate += coeffs[-1] * Basic.log(irreducible)

        h = together(ff - derivation(candidate) / denom)

        numerator = h.as_numer_denom()[0].expand()

        if not isinstance(numerator, Basic.Add):
            numerator = [numerator]

        collected = {}

        for term in numerator:
            coeff, depend = term.as_independent(*V)

            if depend in collected:
                collected[depend] += coeff
            else:
                collected[depend] = coeff

        solutions = solve(collected.values(), coeffs)

        if solutions is None:
            if linear:
                break
            else:
                linear = True
        else:
            break

    if solutions is not None:
        antideriv = candidate.subs_dict(solutions)

        for C in coeffs:
            if C not in solutions:
                antideriv = antideriv.subs(C, S.Zero)

        antideriv = simplify(antideriv.subs_dict(out_terms)).expand()

        if isinstance(antideriv, Basic.Add):
            return Basic.Add(*antideriv.as_coeff_factors()[1])
        else:
            return antideriv
    else:
        if not rewrite:
            return risch_norman(f, x, rewrite=True)
        else:
            return None
