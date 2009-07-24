
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.symbol import Symbol
from sympy.core.basic import Basic, S
from sympy.core.numbers import Integer
from sympy.core.sympify import sympify

from sympy.core.numbers import igcd, ilcm

from polynomial import Poly, SymbolsError, MultivariatePolyError

from monomial import monomial_cmp, monomial_lcm, \
    monomial_gcd, monomial_mul, monomial_div

from sympy.utilities.iterables import all, any

def poly_div(f, g, *symbols):
    """Generalized polynomial division with remainder.

       Given polynomial f and a set of polynomials g = (g_1, ..., g_n)
       compute a set of quotients q = (q_1, ..., q_n) and remainder r
       such that f = q_1*f_1 + ... + q_n*f_n + r, where r = 0 or r is
       a completely reduced polynomial with respect to g.

       In particular g can be a tuple, list or a singleton. All g_i
       and f can be given as Poly class instances or as expressions.

       For more information on the implemented algorithm refer to:

       [1] D. Cox, J. Little, D. O'Shea, Ideals, Varieties and
           Algorithms, Springer, Second Edition, 1997, pp. 62

       [2] I.A. Ajwa, Z. Liu, P.S. Wang, Groebner Bases Algorithm,
           http://citeseer.ist.psu.edu/ajwa95grbner.html, 1995

    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise SymbolsError("Redundant symbols were given")

    f, g = f.unify_with(g)

    symbols, flags = f.symbols, f.flags

    r = Poly((), *symbols, **flags)

    if isinstance(g, Basic):
        if g.is_constant:
            if g.is_zero:
                raise ZeroDivisionError
            elif g.is_one:
                return f, r
            else:
                return f.div_term(g.LC), r

        if g.is_monomial:
            LC, LM = g.lead_term

            q_coeffs, q_monoms = [], []
            r_coeffs, r_monoms = [], []

            for coeff, monom in f.iter_terms():
                quotient = monomial_div(monom, LM)

                if quotient is not None:
                    coeff /= LC

                    q_coeffs.append(Poly.cancel(coeff))
                    q_monoms.append(quotient)
                else:
                    r_coeffs.append(coeff)
                    r_monoms.append(monom)

            return (Poly((q_coeffs, q_monoms), *symbols, **flags),
                    Poly((r_coeffs, r_monoms), *symbols, **flags))

        g, q = [g], [r]
    else:
        q = [r] * len(g)

    while not f.is_zero:
        for i, h in enumerate(g):
            monom = monomial_div(f.LM, h.LM)

            if monom is not None:
                coeff = Poly.cancel(f.LC / h.LC)

                q[i] = q[i].add_term(coeff, monom)
                f -= h.mul_term(coeff, monom)

                break
        else:
            r = r.add_term(*f.LT)
            f = f.kill_lead_term()

    if len(q) != 1:
        return q, r
    else:
        return q[0], r

def poly_quo(f, g, *symbols):
    """Returns polynomial quotient. """
    return poly_div(f, g, *symbols)[0]

def poly_rem(f, g, *symbols):
    """Returns polynomial remainder. """
    return poly_div(f, g, *symbols)[1]

def poly_pdiv(f, g, *symbols):
    """Univariate polynomial pseudo-division with remainder.

       Given univariate polynomials f and g over an integral domain D[x]
       applying classical division algorithm to LC(g)**(d + 1) * f and g
       where  d = max(-1, deg(f) - deg(g)),  compute polynomials q and r
       such that LC(g)**(d + 1)*f = g*q + r and r = 0 or deg(r) < deg(g).
       Polynomials q and r are called the pseudo-quotient of f by g and
       the pseudo-remainder of f modulo g respectively.

       For more information on the implemented algorithm refer to:

       [1] M. Bronstein, Symbolic Integration I: Transcendental
           Functions, Second Edition, Springer-Verlang, 2005

    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise SymbolsError("Redundant symbols were given")

    f, g = f.unify_with(g)

    if f.is_multivariate:
        raise MultivariatePolyError(f)

    symbols, flags = f.symbols, f.flags

    q, r = Poly((), *symbols, **flags), f
    coeff, N = g.LC, f.degree - g.degree + 1

    while not r.is_zero:
        M = r.degree - g.degree

        if M < 0:
            break
        else:
            T, N = (r.LC, (M,)), N - 1

            q = q.mul_term(coeff).add_term(*T)
            r = r.mul_term(coeff)-g.mul_term(*T)

    return (q.mul_term(coeff**N), r.mul_term(coeff**N))

def poly_pquo(f, g, *symbols):
    """Returns polynomial pseudo-quotient. """
    return poly_pdiv(f, g, *symbols)[0]

def poly_prem(f, g, *symbols):
    """Returns polynomial pseudo-remainder. """
    return poly_pdiv(f, g, *symbols)[1]

def poly_groebner(f, *symbols, **flags):
    """Computes reduced Groebner basis for a set of polynomials.

       Given a set of multivariate polynomials F, find another set G,
       such that Ideal F = Ideal G and G is a reduced Groebner basis.

       The resulting basis is unique and has monic generators.

       Groebner bases can be used to choose specific generators for a
       polynomial ideal. Because these bases are unique you can check
       for ideal equality by comparing the Groebner bases.  To see if
       one polynomial lies in an ideal, divide by the elements in the
       base and see if the remainder vanishes.

       They can also be used to  solve systems of polynomial equations
       as,  by choosing lexicographic ordering,  you can eliminate one
       variable at a time, provided that the ideal is zero-dimensional
       (finite number of solutions).

       >>> from sympy import *
       >>> x,y = symbols('xy')

       >>> G = poly_groebner([x**2 + y**3, y**2-x], x, y, order='lex')

       >>> [ g.as_basic() for g in G ]
       [x - y**2, y**3 + y**4]

       For more information on the implemented algorithm refer to:

       [1] N.K. Bose, B. Buchberger, J.P. Guiver, Multidimensional
           Systems Theory and Applications, Springer, 2003, pp. 98+

       [2] A. Giovini, T. Mora, "One sugar cube, please" or Selection
           strategies in Buchberger algorithm, Proc. ISSAC '91, ACM

       [3] I.A. Ajwa, Z. Liu, P.S. Wang, Groebner Bases Algorithm,
           http://citeseer.ist.psu.edu/ajwa95grbner.html, 1995

       [4] D. Cox, J. Little, D. O'Shea, Ideals, Varieties and
           Algorithms, Springer, Second Edition, 1997, pp. 62

    """
    if isinstance(f, (tuple, list, set)):
        f, g = f[0], list(f[1:])

        if not isinstance(f, Poly):
            f = Poly(f, *symbols, **flags)
        elif symbols or flags:
            raise SymbolsError("Redundant symbols or flags were given")

        f, g = f.unify_with(g)

        symbols, flags = f.symbols, f.flags
    else:
        if not isinstance(f, Poly):
            f = Poly(f, *symbols, **flags)
        elif symbols or flags:
            raise SymbolsError("Redundant symbols or flags were given")

        return [f.as_monic()]

    compare = monomial_cmp(flags.get('order'))

    f = [ h for h in [f] + g if h ]

    if not f:
        return [Poly((), *symbols, **flags)]

    R, P, G, B, F = set(), set(), set(), {}, {}

    for i, h in enumerate(f):
        F[h] = i; R.add(i)

    def normal(g, H):
        h = poly_div(g, [ f[i] for i in H ])[1]

        if h.is_zero:
            return None
        else:
            if not h in F:
                F[h] = len(f)
                f.append(h)

            return F[h], h.LM

    def generate(R, P, G, B):
        while R:
            h = normal(f[R.pop()], G | P)

            if h is not None:
                k, LM = h

                G0 = set(g for g in G if monomial_div(f[g].LM, LM))
                P0 = set(p for p in P if monomial_div(f[p].LM, LM))

                G, P, R = G - G0, P - P0 | set([k]), R | G0 | P0

                for i, j in set(B):
                    if i in G0 or j in G0:
                        del B[(i, j)]

        G |= P

        for i in G:
            for j in P:
                if i == j:
                    continue

                if i < j:
                   k = (i, j)
                else:
                   k = (j, i)

                if k not in B:
                    B[k] = monomial_lcm(f[i].LM, f[j].LM)

        G = set([ normal(f[g], G - set([g]))[0] for g in G ])

        return R, P, G, B

    R, P, G, B = generate(R, P, G, B)

    while B:
        k, M = B.items()[0]

        for l, N in B.iteritems():
            if compare(M, N) == 1:
                k, M = l, N

        del B[k]

        i, j = k[0], k[1]
        p, q = f[i], f[j]

        p_LM, q_LM = p.LM, q.LM

        if M == monomial_mul(p_LM, q_LM):
            continue

        criterion = False

        for g in G:
            if g == i or g == j:
                continue

            if (min(i, g), max(i, g)) not in B:
                continue

            if (min(j, g), max(j, g)) not in B:
                continue

            if not monomial_div(M, f[g].LM):
                continue

            criterion = True
            break

        if criterion:
            continue

        p = p.mul_term(1/p.LC, monomial_div(M, p_LM))
        q = q.mul_term(1/q.LC, monomial_div(M, q_LM))

        h = normal(p - q, G)

        if h is not None:
            k, LM = h

            G0 = set(g for g in G if monomial_div(f[g].LM, LM))

            R, P, G = G0, set([k]), G - G0

            for i, j in set(B):
                if i in G0 or j in G0:
                    del B[(i, j)]

            R, P, G, B = generate(R, P, G, B)

    G = [ f[g].as_monic() for g in G ]

    G = sorted(G, compare, lambda p: p.LM)

    return list(reversed(G))

def poly_lcm(f, g, *symbols):
    """Computes least common multiple of two polynomials.

       Given two univariate polynomials,  the LCM is computed  via
       f*g = gcd(f, g)*lcm(f, g) formula. In multivariate case, we
       compute the unique generator of the intersection of the two
       ideals, generated by f and g.  This is done by computing  a
       Groebner basis, with respect to any lexicographic ordering,
       of t*f and (1 - t)*g, where t is an unrelated symbol and
       filtering out solution that does not contain t.

       For more information on the implemented algorithm refer to:

       [1] D. Cox, J. Little, D. O'Shea, Ideals, Varieties and
           Algorithms, Springer, Second Edition, 1997, pp. 187

    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise SymbolsError("Redundant symbols were given")

    f, g = f.unify_with(g)

    symbols, flags = f.symbols, f.flags

    if f.is_monomial and g.is_monomial:
        monom = monomial_lcm(f.LM, g.LM)

        fc, gc = f.LC, g.LC

        if fc.is_Rational and gc.is_Rational:
            coeff = Integer(ilcm(fc.p, gc.p))
        else:
            coeff = S.One

        return Poly((coeff, monom), *symbols, **flags)

    fc, f = f.as_primitive()
    gc, g = g.as_primitive()

    lcm = ilcm(int(fc), int(gc))

    if f.is_multivariate:
        t = Symbol('t', dummy=True)
        lex = { 'order' : 'lex' }

        f_monoms = [ (1,) + monom for monom in f.monoms ]

        F = Poly((f.coeffs, f_monoms), t, *symbols, **lex)

        g_monoms = [ (0,) + monom for monom in g.monoms ] + \
                   [ (1,) + monom for monom in g.monoms ]

        g_coeffs = list(g.coeffs) + [ -coeff for coeff in g.coeffs ]
        G = Poly(dict(zip(g_monoms, g_coeffs)), t, *symbols, **lex)

        def independent(h):
            return all(not monom[0] for monom in h.monoms)

        H = [ h for h in poly_groebner((F, G)) if independent(h) ]

        if lcm != 1:
            h_coeffs = [ coeff*lcm for coeff in H[0].coeffs ]
        else:
            h_coeffs = H[0].coeffs

        h_monoms = [ monom[1:] for monom in H[0].monoms ]

        return Poly(dict(zip(h_monoms, h_coeffs)), *symbols, **flags)
    else:
        h = poly_div(f * g, poly_gcd(f, g))[0]

        if lcm != 1:
            return h.mul_term(lcm / h.LC)
        else:
            return h.as_monic()

def poly_gcd(f, g, *symbols):
    """Compute greatest common divisor of two polynomials.

       Given two univariate polynomials, subresultants are used
       to compute the GCD.  In multivariate case Groebner basis
       approach is used together with f*g = gcd(f, g)*lcm(f, g)
       well known formula.

       For more information on the implemented algorithm refer to:

       [1] D. Cox, J. Little, D. O'Shea, Ideals, Varieties and
           Algorithms, Springer, Second Edition, 1997, pp. 187

    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise SymbolsError("Redundant symbols were given")

    f, g = f.unify_with(g)

    symbols, flags = f.symbols, f.flags

    if f.is_zero and g.is_zero:
        return f

    if f.is_constant:
        if f.is_zero:
            cont, g = g.as_primitive()
            return g.mul_term(cont / g.LC)
        if f.is_one:
            return f

    if g.is_constant:
        if g.is_zero:
            cont, f = f.as_primitive()
            return f.mul_term(cont / f.LC)
        if g.is_one:
            return g

    if f.is_monomial and g.is_monomial:
        monom = monomial_gcd(f.LM, g.LM)

        fc, gc = f.LC, g.LC

        if fc.is_Rational and gc.is_Rational:
            coeff = Integer(igcd(fc.p, gc.p))
        else:
            coeff = S.One

        return Poly((coeff, monom), *symbols, **flags)

    cf, f = f.as_primitive()
    cg, g = g.as_primitive()

    gcd = igcd(int(cf), int(cg))

    if f.is_multivariate:
        h = poly_div(f*g, poly_lcm(f, g))[0]
    else:
        h = poly_subresultants(f, g, res=False)[-1]

    if gcd != 1:
        return h.mul_term(gcd / h.LC)
    else:
        return h.as_monic()

def poly_gcdex(f, g, *symbols):
    """Extended Euclidean algorithm.

       Given univariate polynomials f and g over an Euclidean domain,
       computes polynomials s, t and h,  such that h = gcd(f, g) and
       s*f + t*g = h.

       For more information on the implemented algorithm refer to:

       [1] M. Bronstein, Symbolic Integration I: Transcendental
           Functions, Second Edition, Springer-Verlang, 2005

    """
    s, h = poly_half_gcdex(f, g, *symbols)
    return s, poly_div(h - s*f, g)[0], h

def poly_half_gcdex(f, g, *symbols):
    """Half extended Euclidean algorithm.

       Efficiently computes gcd(f, g)  and one of the coefficients
       in extended Euclidean algorithm. Formally, given univariate
       polynomials f and g over an Euclidean domain, computes s
       and h, such that h = gcd(f, g) and s*f = h (mod g).

       For more information on the implemented algorithm refer to:

       [1] M. Bronstein, Symbolic Integration I: Transcendental
           Functions, Second Edition, Springer-Verlang, 2005

    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise SymbolsError("Redundant symbols were given")

    f, g = f.unify_with(g)

    if f.is_multivariate:
        raise MultivariatePolyError(f)

    symbols, flags = f.symbols, f.flags

    a = Poly(S.One, *symbols, **flags)
    b = Poly((), *symbols, **flags)

    while not g.is_zero:
        q, r = poly_div(f, g)

        f, g = g, r
        c = a - q*b
        a, b = b, c

    return a.div_term(f.LC), f.as_monic()

def poly_resultant(f, g, *symbols):
    """Computes resultant of two univariate polynomials.

       Resultants are a classical algebraic tool for determining if
       a  system of n polynomials in n-1 variables have common root
       without explicitly solving for the roots.

       They are efficiently represented as  determinants of Bezout
       matrices whose entries are computed using O(n**2) additions
       and multiplications where n = max(deg(f), deg(g)).

       >>> from sympy import *
       >>> x,y = symbols('xy')

       Polynomials x**2-1 and (x-1)**2 have common root:

       >>> poly_resultant(x**2-1, (x-1)**2, x)
       0

       For more information on the implemented algorithm refer to:

       [1] Eng-Wee Chionh, Fast Computation of the Bezout and Dixon
           Resultant Matrices, Journal of Symbolic Computation, ACM,
           Volume 33, Issue 1, January 2002, Pages 13-29

    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise SymbolsError("Redundant symbols were given")

    f, g = f.unify_with(g)

    if f.is_multivariate:
        raise MultivariatePolyError(f)

    n, m = f.degree, g.degree

    N = max(n, m)

    if n < m:
        p = f.as_uv_dict()
        q = g.as_uv_dict()
    else:
        q = f.as_uv_dict()
        p = g.as_uv_dict()

    import sympy.matrices

    B = sympy.matrices.zeros(N)

    for i in xrange(N):
        for j in xrange(i, N):
            if i in p and j+1 in q:
                B[i, j] += p[i] * q[j+1]

            if j+1 in p and i in q:
                B[i, j] -= p[j+1] * q[i]

    for i in xrange(1, N-1):
        for j in xrange(i, N-1):
            B[i, j] += B[i-1, j+1]

    for i in xrange(N):
        for j in xrange(i+1, N):
            B[j, i] = B[i, j]

    det = B.det()

    if not det:
        return det
    else:
        if n >= m:
            det /= f.LC**(n-m)
        else:
            det /= g.LC**(m-n)

        sign = (-1)**(n*(n-1)//2)

        if det.is_Atom:
            return sign * det
        else:
            return sign * Poly.cancel(det)

def poly_subresultants(f, g, *symbols, **flags):
    """Computes subresultant PRS of two univariate polynomials.

       Polynomial remainder sequence (PRS) is a fundamental tool in
       computer algebra as it gives as a sub-product the polynomial
       greatest common divisor (GCD), provided that the coefficient
       domain is an unique factorization domain.

       There are several methods for computing PRS, eg.: Euclidean
       PRS, where the most famous algorithm is used, primitive PRS
       and, finally, subresultants which are implemented here.

       The Euclidean approach is reasonably efficient but suffers
       severely from coefficient growth.  The primitive algorithm
       avoids this but requires a lot of coefficient computations.

       Subresultants solve both problems and so it is efficient and
       have moderate coefficient growth. The current implementation
       uses pseudo-divisions  which is well suited for coefficients
       in integral domains or number fields.

       Formally,  given univariate polynomials f and g over an UFD,
       then a sequence (R_0, R_1, ..., R_k, 0, ...) is a polynomial
       remainder sequence where R_0 = f, R_1 = g, R_k != 0 and R_k
       is similar to gcd(f, g).

       The result is returned as tuple (res, R) where R is the PRS
       sequence and res is the resultant of the input polynomials.

       If only polynomial remainder sequence is important,  then by
       setting res=False in keyword arguments expensive computation
       of the resultant can be avoided (only PRS is returned).

       For more information on the implemented algorithm refer to:

       [1] M. Bronstein, Symbolic Integration I: Transcendental
           Functions, Second Edition, Springer-Verlang, 2005

       [2] M. Keber, Division-Free computation of subresultants
           using Bezout matrices, Tech. Report MPI-I-2006-1-006,
           Saarbrucken, 2006

    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise SymbolsError("Redundant symbols were given")

    f, g = f.unify_with(g)

    if f.is_multivariate:
        raise MultivariatePolyError(f)
    else:
        symbols = f.symbols

    n, m = f.degree, g.degree

    if n < m:
        f, g = g, f
        n, m = m, n

    R = [f, g]

    d = n - m

    b = S(-1)**(d + 1)
    c = S(-1)

    B, D = [b], [d]

    h = poly_prem(f, g)
    h = h.mul_term(b)

    while not h.is_zero:
        k = h.degree
        R.append(h)

        lc = g.LC

        C = (-lc)**d / c**(d-1)
        c = Poly.cancel(C)

        b = -lc * c**(m-k)

        f, g, m, d = g, h, k, m-k

        B.append(b)
        D.append(d)

        h = poly_prem(f, g)
        h = h.div_term(b)

    if not flags.get('res', True):
        return R

    if R[-1].degree > 0:
        return (Poly((), *symbols), R)
    if R[-2].is_one:
        return (R[-1], R)

    s, c, i = 1, S(1), 1

    for b, d in zip(B, D)[:-1]:
        u = R[i-1].degree
        v = R[i  ].degree
        w = R[i+1].degree

        if u % 2 and v % 2:
            s = -s

        lc = R[i].LC

        C = c*(b/lc**(1 + d))**v * lc**(u - w)
        c = Poly.cancel(C)

        i += 1

    j = R[-2].degree

    return (R[-1]**j*s*c, R)

def poly_sqf(f, *symbols):
    """Compute square-free decomposition of an univariate polynomial.

       Given an univariate polynomial f over an unique factorization domain
       returns tuple (f_1, f_2, ..., f_n),  where all  A_i are co-prime and
       square-free polynomials and f = f_1 * f_2**2 * ... * f_n**n.

       >>> from sympy import *
       >>> x,y = symbols('xy')

       >>> p, q = poly_sqf(x*(x+1)**2, x)

       >>> p.as_basic()
       x
       >>> q.as_basic()
       1 + x

       For more information on the implemented algorithm refer to:

       [1] M. Bronstein, Symbolic Integration I: Transcendental
           Functions, Second Edition, Springer-Verlang, 2005

       [2] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           Second Edition, Cambridge University Press, 2003

    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise SymbolsError("Redundant symbols were given")

    if f.is_multivariate:
        raise MultivariatePolyError(f)

    coeff, f = f.as_primitive()

    sqf = []

    h = f.diff()

    g = poly_gcd(f, h)

    p = poly_div(f, g)[0]
    q = poly_div(h, g)[0]

    p, q = poly_reduce(p, q)

    while True:
        h = q - p.diff()

        if h.is_zero:
            break

        g = poly_gcd(p, h)

        sqf.append(g)

        p = poly_div(p, g)[0]
        q = poly_div(h, g)[0]

        p, q = poly_reduce(p, q)

    sqf.append(p)

    head, tail = sqf[0], sqf[1:]
    head = head.mul_term(coeff)

    return [head] + tail

def poly_decompose(f, *symbols):
    """Computes functional decomposition of an univariate polynomial.

       Besides factorization and square-free decomposition, functional
       decomposition is another important, but very different,  way of
       breaking down polynomials into simpler parts.

       Formally given an univariate polynomial f with coefficients in a
       field of characteristic zero, returns tuple (f_1, f_2, ..., f_n)
       where f = f_1 o f_2 o ... f_n = f_1(f_2(... f_n)) and f_2, ...,
       f_n are monic and homogeneous polynomials of degree at least 2.

       Unlike factorization, complete functional decompositions of
       polynomials are not unique, consider examples:

        [1] f o g = f(x + b) o (g - b)
        [2] x**n o x**m = x**m o x**n
        [3] T_n o T_m = T_m o T_n

       where T_n and T_m are Chebyshev polynomials.

       >>> from sympy import *
       >>> x,y = symbols('xy')

       >>> p, q = poly_decompose(x**4+2*x**2 + y, x)

       >>> p.as_basic()
       y + 2*x + x**2
       >>> q.as_basic()
       x**2

       For more information on the implemented algorithm refer to:

       [1] D. Kozen, S. Landau, Polynomial decomposition algorithms,
           Journal of Symbolic Computation 7 (1989), pp. 445-456

    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise SymbolsError("Redundant symbols were given")

    if f.is_multivariate:
        raise MultivariatePolyError(f)

    symbols = f.symbols
    flags = f.flags

    def right_factor(f, s):
        n, lc = f.degree, f.LC

        f = f.as_uv_dict()
        q = { s : S.One }

        r = n // s

        for k in xrange(1, s):
            coeff = S.Zero

            for j in xrange(0, k):
                if not n+j-k in f:
                    continue

                if not s-j in q:
                    continue

                fc, qc = f[n+j-k], q[s-j]

                coeff += (k - r*j)*fc*qc

            if coeff is not S.Zero:
                q[s-k] = coeff / (k*r*lc)

        return Poly(q, *symbols, **flags)

    def left_factor(f, h):
        g, i = {}, 0

        while not f.is_zero:
            q, r = poly_div(f, h)

            if not r.is_constant:
                return None
            else:
                if r.LC is not S.Zero:
                    g[i] = r.LC

                f, i = q, i + 1

        return Poly(g, *symbols, **flags)

    def decompose(f):
        deg = f.degree

        for s in xrange(2, deg):
            if deg % s != 0:
                continue

            h = right_factor(f, s)

            if h is not None:
                g = left_factor(f, h)

                if g is not None:
                    return (g, h)

        return None

    F = []

    while True:
        result = decompose(f)

        if result is not None:
            f, h = result
            F = [h] + F
        else:
            break

    return [f] + F

def poly_reduce(f, g, *symbols):
    """Removes common content from a pair of polynomials.

       >>> from sympy import *
       >>> x = Symbol('x')

       >>> f = Poly(2930944*x**6 + 2198208*x**4 + 549552*x**2 + 45796, x)
       >>> g = Poly(17585664*x**5 + 8792832*x**3 + 1099104*x, x)

       >>> F, G = poly_reduce(f, g)

       >>> F
       Poly(64*x**6 + 48*x**4 + 12*x**2 + 1, x)
       >>> G
       Poly(384*x**5 + 192*x**3 + 24*x, x)

    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise SymbolsError("Redundant symbols were given")

    f, g = f.unify_with(g)

    fc = int(f.content)
    gc = int(g.content)

    cont = igcd(fc, gc)

    if cont != 1:
        f = f.div_term(cont)
        g = g.div_term(cont)

    return f, g

def poly_discriminant(p):
    """
    Returns the discriminant of a polynomial.

    The discriminant of a univariate polynomial p of degree n is defined as
    n**(n*(n-1)/2)/a_n*resultant(p, p'), where p' is the derivative of p and a_n
    is the leading coefficient of p.  Because the resultant of two polynomials
    vanishes identically whenever the two polynomials share a root, and a
    polynomial shares a root with its derivative if and only if the root is a
    repeated root, it follows that the discriminant of a polynomial vanishes
    identically if and only if the polynomial has a repeated root.

    See also:
    <http://en.wikipedia.org/wiki/Discriminant>

    Example:
    >>> from sympy import *
    >>> a, b, c, x = symbols('abcx')
    >>> discriminant(Poly(a*x**2 + b*x + c, x))
    -4*a*c + b**2
    >>> discriminant(Poly(2*x**5 + x**2 + 10, x))
    500004320
    >>> discriminant(Poly((x-1)*(x+1), x))
    4
    >>> discriminant(Poly((x-1)**2*(x+1), x))
    0
    """
    if len(p.symbols) != 1:
        raise NotImplementedError("Only univariate polynomials are supported.")
    def negonetox(x):
        """(-1)**x"""
        return -2*(x%2) + 1
    if p.degree == 0:
        return S(0)
    return (negonetox(S(p.degree)*(S(p.degree) - 1)/2)/p.lead_coeff*\
    poly_resultant(p, p.diff(p.symbols[0]))).expand()

