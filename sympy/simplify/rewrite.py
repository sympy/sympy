"""Module 'rewrite.py' contains advanced term rewriting methods concerning
   partial fraction decomposition, combinig together and collecting terms.
"""

from sympy.core import Basic, S, Symbol, Add, Function
from sympy.core.methods import NoRelMeths, ArithMeths

from sympy.polynomials import div, quo, rem, gcd
from sympy.simplify import normal, together
from sympy.polynomials.factor_ import sqf

###############################################################################
### THIS HAS TO BE MOVED TO APPROPRIATE MODULE

def ext_gcd(p, q, x):
    U = (p, S.One, S.Zero)
    V = (q, S.Zero, S.One)

    while True:
        q = quo(U[0], V[0], x)

        U, V = V, [ (a - q*b).expand() for a, b in zip(U, V) ]

        if isinstance(V[0], Basic.Zero):
            return U

class RootOf(Basic):

    precedence = Basic.Apply_precedence

    def __new__(cls, f, z, **assumptions):
        f = Basic.sympify(f)

        if not f.is_polynomial(z):
            return f

        obj = Basic.__new__(cls, **assumptions)
        obj._args = (f.as_polynomial(z), z)

        return obj

    def tostr(self, level=0):
        return 'RootOf%s' % str( (self._args[0].as_basic(), self._args[1]) )

    def doit(self, **hints):
        return self

        #if not hints.get('roots', True):
        #    return self
        #else:
        #    domain = hints.get('domain', None)
        #
        #    if domain is None:
        #        domain = self.domain
        #
        #    return {} # TBD


###############################################################################

def apart(f, z, domain=None, index=None):
    """Computes full partial fraction decomposition of a univariate
       rational function over the algebraic closure of its field of
       definition. Although only gcd operations over the initial
       field are required, the expansion is returned in a formal
       form with linear denominators.

       However it is possible to force expansion of the resulting
       formal summations, and so factorization over a specified
       domain is performed.

       To specify the desired behavior of the algorithm use the
       'domain' keyword. Setting it to None, which is done be
       default, will result in no factorization at all.

       Otherwise it can be assigned with one of Z, Q, R, C domain
       specifiers and the formal partial fraction expansion will
       be rewritten using all possible roots over this domain.

       If the resulting expansion contains formal summations, then
       for all those a single dummy index variable named 'a' will
       be generated. To change this default behavior issue new
       name via 'index' keyword.

       For more information on the implemented algorithm refer to:

       [1] M. Bronstein, B. Salvy, Full partial fraction decomposition
           of rational functions, in: M. Bronstein, ed., Proceedings
           ISSAC '93, ACM Press, Kiev, Ukraine, 1993, pp. 157-160.

    """
    f = Basic.sympify(f)

    if isinstance(f, Basic.Add):
        return Add(*[ apart(g) for g in f ])
    else:
        if f.is_fraction(z):
            f = normal(f, z)
        else:
            return f

        P, Q = f.as_numer_denom()

        if not Q.has(z):
            return f

        u = Function('u')(z)

        if index is None:
            A = Symbol('a', dummy=True)
        else:
            A = Symbol(index)

        partial, r = div(P, Q, z)
        f, q, U = r / Q, Q, []

        for k, d in enumerate(sqf(q, z)):
            n, d = k + 1, d.as_basic()
            U += [ u.diff(z, k) ]

            h = normal(f * d**n, z) / u**n

            H, subs = [h], []

            for j in range(1, n):
                H += [ H[-1].diff(z) / j ]

            for j in range(1, n+1):
                subs += [ (U[j-1], d.diff(z, j) / j) ]

            for j in range(0, n):
                P, Q = together(H[j]).as_numer_denom()

                for i in range(0, j+1):
                    P = P.subs(*subs[j-i])

                Q = Q.subs(*subs[0])

                G = gcd(P, d, z)
                D = quo(d, G, z)

                g, B, _ = ext_gcd(Q, D, z)
                b = rem(P * B / g, D, z)

                term = b.subs(z, A) / (z - A)**(n-j)

                if domain is None:
                    a = D.diff(z)

                    if not a.has(z):
                        partial += term.subs(A, -D.subs(z, 0) / a)
                    else:
                        partial += Basic.Sum(term, (A, Basic.RootOf(D, z)))
                else:
                    raise NotImplementedError

        return partial
