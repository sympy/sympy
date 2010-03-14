"""Partial fraction decomposition algorithms. """

from sympy.core import sympify, Symbol, Function, Lambda
from sympy.polys import Poly, RootSum, div, cancel
from sympy.utilities import threaded

@threaded()
def apart(f, z=None, **args):
    """Computes partial fraction decomposition of a rational function.

       Given a rational function 'f', performing only gcd operations
       over the algebraic closure of the initial field of definition,
       compute full partial fraction decomposition with fractions
       having linear denominators.

       For all other kinds of expressions the input is returned in an
       unchanged form. Note however, that `apart` function can thread
       over sums and relational operators.

       Note that no factorization of the initial denominator of `f` is
       needed.  The final decomposition is formed in terms of a sum of
       RootSum instances.  By default RootSum tries to compute all its
       roots to simplify itself. This behavior can be however avoided
       by setting the keyword flag evaluate=False, which will make this
       function return a formal decomposition.

           >>> from sympy import apart
           >>> from sympy.abc import x, y

           >>> apart(y/(x+2)/(x+1), x)
           y/(1 + x) - y/(2 + x)

           >>> apart(1/(1+x**5), x, evaluate=False)
           RootSum(Lambda(_a, -_a/(5*(x - _a))), x**5 + 1, x, domain='ZZ')

       References
       ==========

       .. [Bronstein93] M. Bronstein, B. Salvy, Full partial fraction
           decomposition of rational functions, Proceedings ISSAC '93,
           ACM Press, Kiev, Ukraine, 1993, pp. 157-160.

    """
    f = cancel(f)

    if z is None:
        symbols = f.atoms(Symbol)

        if not symbols:
            return f

        if len(symbols) == 1:
            z = list(symbols)[0]
        else:
            raise ValueError("multivariate partial fractions are not supported")

    P, Q = f.as_numer_denom()

    if not Q.has(z):
        return f

    partial, r = div(P, Q, z)
    f, q, U = r / Q, Q, []

    u = Function('u')(z)
    a = Symbol('a', dummy=True)

    q = Poly(q, z)

    for d, n in q.sqf_list(all=True, include=True):
        b = d.as_basic()
        U += [ u.diff(z, n-1) ]

        h = cancel(f*b**n) / u**n

        H, subs = [h], []

        for j in range(1, n):
            H += [ H[-1].diff(z) / j ]

        for j in range(1, n+1):
            subs += [ (U[j-1], b.diff(z, j) / j) ]

        for j in range(0, n):
            P, Q = cancel(H[j]).as_numer_denom()

            for i in range(0, j+1):
                P = P.subs(*subs[j-i])

            Q = Q.subs(*subs[0])

            P = Poly(P, z)
            Q = Poly(Q, z)

            G = P.gcd(d)
            D = d.exquo(G)

            B, g = Q.half_gcdex(D)
            b = (P * B.exquo(g)).rem(D)

            numer = b.as_basic()
            denom = (z-a)**(n-j)

            expr = numer.subs(z, a) / denom

            partial += RootSum(Lambda(a, expr), D, **args)

    return partial

