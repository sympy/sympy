
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.symbol import Symbol
from sympy.core.numbers import Integer
from sympy.core.sympify import sympify
from sympy.core.basic import Basic, S, C, Atom

from polynomial import Poly, PolynomialError
from monomial import monomial_lcm, monomial_div

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

    if isinstance(g, (tuple, list)):
        q = [r] * len(g)
    else:
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

def poly_groebner(f, *symbols):
    pass # TBD : groebner needs sugar

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
