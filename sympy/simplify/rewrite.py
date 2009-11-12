"""Module 'rewrite.py' contains advanced term rewriting methods concerning
   partial fraction decomposition, combinig together and collecting terms.
"""

from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.relational import Relational
from sympy.core.symbol import Symbol, Temporary
from sympy.core.numbers import Integer, Rational
from sympy.core.function import Function, Lambda
from sympy.core.basic import Basic, S, C, Atom, sympify

from sympy.utilities import threaded
from sympy.simplify import together
from sympy.functions import exp

from sympy.polys import Poly, RootSum, div, quo, gcd
from sympy.polys.algorithms import poly_quo, poly_rem, \
    poly_sqf, poly_gcd, poly_half_gcdex

@threaded()
def cancel(f, *symbols):
    """Cancel common factors in a given rational function.

       Given a quotient of polynomials, performing only gcd and quo
       operations in polynomial algebra,  return rational function
       with numerator and denominator of minimal total degree in
       an expanded form.

       For all other kinds of expressions the input is returned in
       an unchanged form. Note however, that 'cancel' function can
       thread over sums and relational operators.

       Additionally you can specify a list of variables to perform
       cancelation more efficiently using only those symbols.

       >>> from sympy import cancel, sqrt
       >>> from sympy.abc import x, y

       >>> cancel((x**2-1)/(x-1))
       1 + x

       >>> cancel((x**2-y**2)/(x-y), x)
       x + y

       >>> cancel((x**2-2)/(x+sqrt(2)))
       x - 2**(1/2)

    """
    return Poly.cancel(f, *symbols)

def trim(f, *symbols, **flags):
    """Cancel common factors in a given formal rational expression.

       Given an arbitrary expression, map all functional components
       to temporary symbols, rewriting this expression to rational
       function form and perform cancelation of common factors.

       When given a rational function or a list of symbols discards
       all functional components, then this procedure is equivalent
       to cancel().

       Note that this procedure can thread over composite objects
       like big operators, matrices, relational operators etc. It
       can be also called recursively (to change this behaviour
       unset 'recursive' flag).

       >>> from sympy import Function, trim, sin

       >>> from sympy.abc import x, y
       >>> f = Function('f')

       >>> trim((f(x)**2+f(x))/f(x))
       1 + f(x)

       >>> trim((x**2+x)/x)
       1 + x

       Recursively simplify expressions:

       >>> trim(sin((f(x)**2+f(x))/f(x)))
       sin(1 + f(x))

    """
    f = sympify(f)

    if isinstance(f, Relational):
        return Relational(trim(f.lhs, *symbols, **flags),
                          trim(f.rhs, *symbols, **flags), f.rel_op)
    #elif isinstance(f, Matrix):
    #    return f.applyfunc(lambda g: trim(g, *symbols, **flags))
    else:
        recursive = flags.get('recursive', True)

        def is_functional(g):
            return not (g.is_Atom or g.is_number) \
                and (not symbols or g.has(*symbols))

        def components(g):
            result = set()

            if is_functional(g):
                if g.is_Add or g.is_Mul:
                    args = []

                    for h in g.args:
                        h, terms = components(h)

                        result |= terms
                        args.append(h)

                    g = g.__class__(*args)
                elif g.is_Pow:
                    if recursive:
                        base = trim(g.base, *symbols, **flags)
                    else:
                        base = g.base

                    if g.exp.is_Rational:
                        if g.exp.is_Integer:
                            if g.exp is S.NegativeOne:
                                h, terms = components(base)
                                return h**S.NegativeOne, terms
                            else:
                                h = base
                        else:
                            h = base**Rational(1, g.exp.q)

                        g = base**g.exp
                    else:
                        if recursive:
                            h = g = base**trim(g.exp, *symbols, **flags)
                        else:
                            h = g = base**g.exp

                    if is_functional(h):
                        result.add(h)
                else:
                    if not recursive:
                        result.add(g)
                    else:
                        g = g.__class__(*[trim(h, *symbols, **flags)
                                          for h in g.args])

                        if is_functional(g):
                            result.add(g)

            return g, result

        if not isinstance(f, Basic) \
           or f.is_number \
           or not f.has_any_symbols(*symbols):
            return f

        f = together(f.expand())
        f, terms = components(f)

        if not terms:
            return Poly.cancel(f, *symbols)
        else:
            mapping, reverse = {}, {}

            for g in terms:
                mapping[g] = Temporary()
                reverse[mapping[g]] = g

            p, q = f.as_numer_denom()
            f = p.expand()/q.expand()

            if not symbols:
                symbols = tuple(f.atoms(Symbol))

            symbols = tuple(mapping.values()) + symbols

            H = Poly.cancel(f.subs(mapping), *symbols)

            if not flags.get('extract', True):
                return H.subs(reverse)
            else:
                def extract(f):
                    p = f.args[0]

                    for q in f.args[1:]:
                        p = gcd(p, q, *symbols)

                        if p.is_number:
                            return S.One, f

                    return p, Add(*[quo(g, p, *symbols) for g in f.args])

                P, Q = H.as_numer_denom()

                if P.is_Add:
                    GP, P = extract(P)
                else:
                    GP = S.One

                if Q.is_Add:
                    GQ, Q = extract(Q)
                else:
                    GQ = S.One

                return ((GP*P)/(GQ*Q)).subs(reverse)

@threaded()
def apart(f, z, **flags):
    """Compute partial fraction decomposition of a rational function.

       Given a rational function 'f', performing only gcd operations
       over the algebraic closue of the initial field of definition,
       compute full partial fraction decomposition with fractions
       having linear denominators.

       For all other kinds of expressions the input is returned in an
       unchanged form. Note however, that 'apart' function can thread
       over sums and relational operators.

       Note that no factorization of the initial denominator of 'f' is
       needed.  The final decomposition is formed in terms of a sum of
       RootSum instances.  By default RootSum tries to compute all its
       roots to simplify itself. This behaviour can be however avoided
       by seting the keyword flag evaluate=False, which will make this
       function return a formal decomposition.

       >>> from sympy import apart
       >>> from sympy.abc import x, y

       >>> apart(y/(x+2)/(x+1), x)
       y/(1 + x) - y/(2 + x)

       >>> apart(1/(1+x**5), x, evaluate=False)
       RootSum(Lambda(_a, -_a/(5*(x - _a))), x**5 + 1, x)

       For more information on the implemented algorithm refer to:

       [1] M. Bronstein, B. Salvy, Full partial fraction decomposition
           of rational functions,  in: M. Bronstein,  ed., Proceedings
           ISSAC '93, ACM Press, Kiev, Ukraine, 1993, pp. 157-160.

    """
    if not f.has(z):
        return f

    f = Poly.cancel(f, z)

    P, Q = f.as_numer_denom()

    if not Q.has(z):
        return f

    partial, r = div(P, Q, z)
    f, q, U = r / Q, Q, []

    u = Function('u')(z)
    a = Symbol('a', dummy=True)

    for k, d in enumerate(poly_sqf(q, z)):
        n, b = k + 1, d.as_basic()
        U += [ u.diff(z, k) ]

        h = together(Poly.cancel(f*b**n, z) / u**n)

        H, subs = [h], []

        for j in range(1, n):
            H += [ H[-1].diff(z) / j ]

        for j in range(1, n+1):
            subs += [ (U[j-1], b.diff(z, j) / j) ]

        for j in range(0, n):
            P, Q = together(H[j]).as_numer_denom()

            for i in range(0, j+1):
                P = P.subs(*subs[j-i])

            Q = Q.subs(*subs[0])

            P, Q = Poly(P, z), Poly(Q, z)

            G = poly_gcd(P, d)
            D = poly_quo(d, G)

            B, g = poly_half_gcdex(Q, D)
            b = poly_rem(P * poly_quo(B, g), D)

            numer = b.as_basic()
            denom = (z-a)**(n-j)

            expr = numer.subs(z, a) / denom

            partial += RootSum(Lambda(a, expr), D, **flags)

    return partial
