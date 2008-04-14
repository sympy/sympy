
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.symbol import Symbol
from sympy.core.numbers import Integer
from sympy.core.sympify import sympify
from sympy.core.basic import Basic, S, C, Atom

from polynomial import Poly, PolynomialError
from monomial import monomial_cmp, monomial_lcm, \
                     monomial_mul, monomial_div

from sympy.simplify import cancel # TBD : move cancel() here
from sympy.matrices import zero

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
        raise PolynomialError

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
                LC = g.lead_coeff

                coeffs = [ coeff / LC for coeff in f.coeffs ]

                for i, coeff in enumerate(coeffs):
                    if not coeff.is_Atom:
                        coeffs[i] = cancel(coeff)

                return Poly((coeffs, f.monoms), *symbols, **flags), r

        if g.is_monomial:
            LC, LM = g.lead_term

            q_coeffs, q_monoms = [], []
            r_coeffs, r_monoms = [], []

            for coeff, monom in f.iter_terms():
                M = monomial_div(monom, LM)

                if M is not None:
                    coeff /= LC

                    if coeff.is_Atom:
                        q_coeffs.append(coeff)
                    else:
                        q_coeffs.append(cancel(coeff))

                    q_monoms.append(M)
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
            M = monomial_div(f.LM, h.LM)

            if M is not None:
                coeff = f.LC / h.LC

                if coeff.is_Atom:
                    T = coeff, M
                else:
                    T = cancel(coeff), M

                P = Poly(T, *symbols, **flags)

                q[i] = q[i].add_term(*T)
                f -= h * P

                break
        else:
            r = r.add_term(*f.LT)
            f = f.kill_lead_term()

    if len(q) != 1:
        return tuple(q), r
    else:
        return q[0], r

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
        raise PolynomialError

    f, g = f.unify_with(g)

    if f.is_multivariate:
        raise PolynomialError

    symbols, flags = f.symbols, f.flags

    q, r = Poly((), *symbols, **flags), f
    coeff, N = g.LC, f.degree - g.degree + 1

    while not r.is_zero:
        M = r.degree - g.degree

        if M < 0:
            break
        else:
            T = Poly((r.LC, (M,)), *symbols, **flags)
            q, r, N = q*coeff + T, r*coeff - g*T, N-1

    return (q * coeff**N, r * coeff**N)

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
            raise PolynomialError

        f, g = f.unify_with(g)

        symbols, flags = f.symbols, f.flags
    else:
        if not isinstance(f, Poly):
            f = Poly(f, *symbols, **flags)
        elif symbols or flags:
            raise PolynomialError

        return (f.as_monic(),)

    compare = monomial_cmp(flags.get('order'))

    f = [ h for h in [f] + g if h ]

    if not f:
        return (Poly((), *symbols, **flags),)

    R, P, G, B, F = set(), set(), set(), {}, {}

    for i, h in enumerate(f):
        F[h] = i; R.add(i)

    def normal(g, H):
        h = poly_div(g, [ f[i] for i in H ])[1]

        if h.is_zero:
            return None
        else:
            if not F.has_key(h):
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

                if not B.has_key(k):
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

            if not B.has_key((min(i, g), max(i, g))):
                continue

            if not B.has_key((min(j, g), max(j, g))):
                continue

            if not monomial_div(M, f[g].LM):
                continue

            criterion = True
            break

        if criterion:
            continue

        p_M = monomial_div(M, p_LM)
        q_M = monomial_div(M, q_LM)

        p *= Poly(((1/p.LC,), (p_M,)), *symbols, **flags)
        q *= Poly(((1/q.LC,), (q_M,)), *symbols, **flags)

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

    return tuple(reversed(G))

def poly_lcm(f, g, *symbols):
    pass # TBD : lcm needs groebner

def poly_gcd(f, g, *symbols):
    pass # TBD : gcd needs groebner

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
        raise PolynomialError

    f, g = f.unify_with(g)

    if f.is_multivariate:
        raise PolynomialError

    symbols, flags = f.symbols, f.flags

    a = Poly(S.One, *symbols, **flags)
    b = Poly((), *symbols, **flags)

    while not g.is_zero:
        q, r = poly_div(f, g)

        f, g = g, r
        c = a - q*b
        a, b = b, c

    return a, f

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
        raise PolynomialError

    f, g = f.unify_with(g)

    if f.is_multivariate:
        raise PolynomialError

    n, m = f.degree, g.degree

    N = max(n, m)

    if n < m:
        p = f.as_uv_dict()
        q = g.as_uv_dict()
    else:
        q = f.as_uv_dict()
        p = g.as_uv_dict()

    B = zero(N)

    for i in xrange(N):
        for j in xrange(i, N):
            if p.has_key(i) and q.has_key(j+1):
                B[i, j] += p[i] * q[j+1]

            if p.has_key(j+1) and q.has_key(i):
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

        sign = (-1)**(n*(n-1)/2)

        if det.is_Atom:
            return sign * det
        else:
            return sign * cancel(det)
