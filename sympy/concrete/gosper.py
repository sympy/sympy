from sympy.core.basic import Basic, S
from sympy.core.symbol import Symbol
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core import sympify

from sympy.polys import gcd, quo, roots, resultant

def normal(f, g, n=None):
    """Given relatively prime univariate polynomials 'f' and 'g',
       rewrite their quotient to a normal form defined as follows:

                       f(n)       A(n) C(n+1)
                       ----  =  Z -----------
                       g(n)       B(n)  C(n)

       where Z is arbitrary constant and A, B, C are monic
       polynomials in 'n' with follwing properties:

           (1) gcd(A(n), B(n+h)) = 1 for all 'h' in N
           (2) gcd(B(n), C(n+1)) = 1
           (3) gcd(A(n), C(n)) = 1

       This normal form, or rational factorization in other words,
       is crucial step in Gosper's algorithm and in difference
       equations solving. It can be also used to decide if two
       hypergeometric are similar or not.

       This procedure will return return triple containig elements
       of this factorization in the form (Z*A, B, C). For example:

       >>> from sympy import Symbol, normal
       >>> n = Symbol('n', integer=True)

       >>> normal(4*n+5, 2*(4*n+1)*(2*n+3), n)
       (1/4, 3/2 + n, 1/4 + n)

    """
    f, g = map(sympify, (f, g))

    p = f.as_poly(n)
    q = g.as_poly(n)

    a, p = p.LC, p.as_monic()
    b, q = q.LC, q.as_monic()

    A = p.as_basic()
    B = q.as_basic()

    C, Z = S.One, a / b

    h = Symbol('h', dummy=True)

    res = resultant(A, B.subs(n, n+h), n)

    nni_roots = roots(res, h, domain='Z',
        predicate=lambda r: r >= 0).keys()

    if not nni_roots:
        return (f, g, S.One)
    else:
        for i in sorted(nni_roots):
            d = gcd(A, B.subs(n, n+i), n)

            A = quo(A, d, n)
            B = quo(B, d.subs(n, n-i), n)

            C *= Mul(*[ d.subs(n, n-j) for j in xrange(1, i+1) ])

        return (Z*A, B, C)

def gosper(term, k, a, n):
    from sympy.solvers import rsolve_poly

    if not term:
        return None
    else:
        p, q = term.as_numer_denom()
        A, B, C = normal(p, q, k)

        B = B.subs(k, k-1)

        R = rsolve_poly([-B, A], C, k)
        symbol = []

        if not (R is None  or  R is S.Zero):
            if symbol != []:
                symbol = symbol[0]

                W = R.subs(symbol, S.Zero)

                if W is S.Zero:
                    R = R.subs(symbol, S.One)
                else:
                    R = W

            Z = B*R*term/C
            return simplify(Z.subs(k, n+1) - Z.subs(k, a))
        else:
            return None
