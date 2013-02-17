"""Algorithms for partial fraction decomposition of rational functions. """

from sympy.polys import Poly, RootSum, cancel, factor
from sympy.polys.polytools import parallel_poly_from_expr
from sympy.polys.polyoptions import allowed_flags, set_defaults

from sympy.core import S, Add, sympify, Function, Lambda, Dummy
from sympy.utilities import numbered_symbols, take, threaded


@threaded
def apart(f, x=None, full=False, **options):
    """
    Compute partial fraction decomposition of a rational function.

    Given a rational function ``f`` compute partial fraction decomposition
    of ``f``. Two algorithms are available: one is based on undetermined
    coefficients method and the other is Bronstein's full partial fraction
    decomposition algorithm.

    Examples
    ========

    >>> from sympy.polys.partfrac import apart
    >>> from sympy.abc import x, y

    By default, using the undetermined coefficients method:
    >>> apart(y/(x + 2)/(x + 1), x)
    -y/(x + 2) + y/(x + 1)

    You can choose the other algorithm by setting full=True:
    >>> apart(y/(x**2 + x + 1), x)
    y/(x**2 + x + 1)
    >>> apart(y/(x**2 + x + 1), x, full=True)
    RootSum(x**2 + x + 1, Lambda(_a, (-2*_a*y/3 - y/3)/(-_a + x)))

    """
    allowed_flags(options, [])

    f = sympify(f)

    if f.is_Atom:
        return f
    else:
        P, Q = f.as_numer_denom()

    options = set_defaults(options, extension=True)
    (P, Q), opt = parallel_poly_from_expr((P, Q), x, **options)

    if P.is_multivariate:
        raise NotImplementedError(
            "multivariate partial fraction decomposition")

    common, P, Q = P.cancel(Q)

    poly, P = P.div(Q, auto=True)
    P, Q = P.rat_clear_denoms(Q)

    if Q.degree() <= 1:
        partial = P/Q
    else:
        if not full:
            partial = apart_undetermined_coeffs(P, Q)
        else:
            partial = apart_full_decomposition(P, Q)

    terms = S.Zero

    for term in Add.make_args(partial):
        if term.has(RootSum):
            terms += term
        else:
            terms += factor(term)

    return common*(poly.as_expr() + terms)


def apart_undetermined_coeffs(P, Q):
    """Partial fractions via method of undetermined coefficients. """
    X = numbered_symbols(cls=Dummy)
    partial, symbols = [], []

    _, factors = Q.factor_list()

    for f, k in factors:
        n, q = f.degree(), Q

        for i in xrange(1, k + 1):
            coeffs, q = take(X, n), q.quo(f)
            partial.append((coeffs, q, f, i))
            symbols.extend(coeffs)

    dom = Q.get_domain().inject(*symbols)
    F = Poly(0, Q.gen, domain=dom)

    for i, (coeffs, q, f, k) in enumerate(partial):
        h = Poly(coeffs, Q.gen, domain=dom)
        partial[i] = (h, f, k)
        q = q.set_domain(dom)
        F += h*q

    system, result = [], S(0)

    for (k,), coeff in F.terms():
        system.append(coeff - P.nth(k))

    from sympy.solvers import solve
    solution = solve(system, symbols)

    for h, f, k in partial:
        h = h.as_expr().subs(solution)
        result += h/f.as_expr()**k

    return result


def apart_full_decomposition(P, Q):
    """
    Bronstein's full partial fraction decomposition algorithm.

    Given a univariate rational function ``f``, performing only GCD
    operations over the algebraic closure of the initial ground domain
    of definition, compute full partial fraction decomposition with
    fractions having linear denominators.

    Note that no factorization of the initial denominator of ``f`` is
    performed. The final decomposition is formed in terms of a sum of
    :class:`RootSum` instances.

    References
    ==========

    1. [Bronstein93]_

    """
    f, x, U = P/Q, P.gen, []

    u = Function('u')(x)
    a = Dummy('a')

    partial = S(0)

    for d, n in Q.sqf_list_include(all=True):
        b = d.as_expr()
        U += [ u.diff(x, n - 1) ]

        h = cancel(f*b**n) / u**n

        H, subs = [h], []

        for j in range(1, n):
            H += [ H[-1].diff(x) / j ]

        for j in range(1, n + 1):
            subs += [ (U[j - 1], b.diff(x, j) / j) ]

        for j in range(0, n):
            P, Q = cancel(H[j]).as_numer_denom()

            for i in range(0, j + 1):
                P = P.subs(*subs[j - i])

            Q = Q.subs(*subs[0])

            P = Poly(P, x)
            Q = Poly(Q, x)

            G = P.gcd(d)
            D = d.quo(G)

            B, g = Q.half_gcdex(D)
            b = (P * B.quo(g)).rem(D)

            numer = b.as_expr()
            denom = (x - a)**(n - j)

            func = Lambda(a, numer.subs(x, a)/denom)
            partial += RootSum(D, func, auto=False)

    return partial


def apart_structured(f, x=None, **options):
    """
    Compute partial fraction decomposition of a rational function
    and return the result in structured form.

    Given a rational function ``f`` compute partial fraction decomposition
    of ``f``. Two algorithms are available: one is based on undetermined
    coefficients method and the other is Bronstein's full partial fraction
    decomposition algorithm.

    Examples
    ========

    >>> from sympy.polys.partfrac import apart_structured, assemble_partfrac_full
    >>> from sympy.abc import x, y

    >>> f = 36 / (x**5 - 2*x**4 - 2*x**3 + 4*x**2 + x - 2)
    >>> pfd = apart_structured(f)
    >>> pfd
    (1,
    Poly(0, x, domain='ZZ'),
    [(Poly(x - 2, x, domain='ZZ'), Lambda(_a, 4), 1),
    (Poly(x**2 - 1, x, domain='ZZ'), Lambda(_a, -3*_a - 6), 2),
    (Poly(x + 1, x, domain='ZZ'), Lambda(_a, -4), 1)])

    >>> assemble_partfrac_full(pfd, x)
    -4/(x + 1) - 3/(x + 1)**2 - 9/(x - 1)**2 + 4/(x - 2)

    See also
    ========

    apart
    """
    allowed_flags(options, [])

    f = sympify(f)

    if f.is_Atom:
        return f
    else:
        P, Q = f.as_numer_denom()

    options = set_defaults(options, extension=True)
    (P, Q), opt = parallel_poly_from_expr((P, Q), x, **options)

    if P.is_multivariate:
        raise NotImplementedError(
            "multivariate partial fraction decomposition")

    common, P, Q = P.cancel(Q)

    poly, P = P.div(Q, auto=True)
    P, Q = P.rat_clear_denoms(Q)

    if Q.degree() <= 1:
        polypart = P/Q
        partial = (common, polypart, ())
    else:
        polypart = poly
        rationalpart = apart_full_decomposition_structured(P, Q)
        partial = (common, polypart, rationalpart)

    return partial


# This could supersede the above full decomposition function
def apart_full_decomposition_structured(P, Q):
    f, x, U = P/Q, P.gen, []

    u = Function('u')(x)
    a = Dummy('a')

    partial = []

    for d, n in Q.sqf_list_include(all=True):
        b = d.as_expr()
        U += [ u.diff(x, n - 1) ]

        h = cancel(f*b**n) / u**n

        H, subs = [h], []

        for j in range(1, n):
            H += [ H[-1].diff(x) / j ]

        for j in range(1, n + 1):
            subs += [ (U[j - 1], b.diff(x, j) / j) ]

        for j in range(0, n):
            P, Q = cancel(H[j]).as_numer_denom()

            for i in range(0, j + 1):
                P = P.subs(*subs[j - i])

            Q = Q.subs(*subs[0])

            P = Poly(P, x)
            Q = Poly(Q, x)

            G = P.gcd(d)
            D = d.quo(G)

            B, g = Q.half_gcdex(D)
            b = (P * B.quo(g)).rem(D)

            numer = b.as_expr()
            denom = (x - a)**(n - j)

            part = (D, Lambda(a, numer.subs(x,a)), n-j)
            partial.append(part)

    return partial


# This could be called internally from apart(...) to get unstructured result
def assemble_partfrac_full(partial_structured, x):
    r"""Reassemble a full partial fraction decomposition
    from a structured result.

    Examples
    ========

    This example is taken from Bronstein's original paper.

    >>> from sympy.polys.partfrac import apart_structured, assemble_partfrac_full
    >>> from sympy.abc import x, y

    >>> f = 36 / (x**5 - 2*x**4 - 2*x**3 + 4*x**2 + x - 2)
    >>> pfd = apart_structured(f)
    >>> pfd
    (1,
    Poly(0, x, domain='ZZ'),
    [(Poly(x - 2, x, domain='ZZ'), Lambda(_a, 4), 1),
    (Poly(x**2 - 1, x, domain='ZZ'), Lambda(_a, -3*_a - 6), 2),
    (Poly(x + 1, x, domain='ZZ'), Lambda(_a, -4), 1)])

    >>> assemble_partfrac_full(pfd, x)
    -4/(x + 1) - 3/(x + 1)**2 - 9/(x - 1)**2 + 4/(x - 2)

    See also
    ========

    apart
    """
    # What is the best way to get x from the data given?

    # Common factor
    common = partial_structured[0]

    # Polynomial part
    polypart = partial_structured[1]
    pfd = polypart.as_expr()

    # Rational parts
    for poly, nf, ex in partial_structured[2]:
        a, nu = nf.args
        func = Lambda(a, nu/(x-a[0])**ex)
        # We need an option to disallow resolution of rootsums even in quadratic case if desired
        pfd += RootSum(poly, func, auto=False)

    return common*pfd
