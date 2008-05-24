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

from sympy.polynomials import factor_, div, quo, rem, gcd, egcd
from sympy.simplify import together
from sympy.matrices import Matrix
from sympy.functions import exp

from sympy.polys import RootSum

def cancel(f, *syms, **flags):
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

       >>> from sympy import *
       >>> x,y = symbols('xy')

       >>> cancel((x**2-1)/(x-1))
       1 + x

       >>> cancel((x**2-y**2)/(x-y), x)
       x + y

       >>> cancel((x**2-2)/(x+sqrt(2)))
       x - 2**(1/2)

    """
    f = sympify(f)

    if syms and not f.has(*syms):
        return f
    elif f.is_Add:
        return Add(*[ cancel(g, *syms, **flags) for g in f.args ])
    elif isinstance(f, Relational):
        return Relational(cancel(f.lhs, *syms, **flags),
            cancel(f.rhs, *syms, **flags), f.rel_op)
    else:
        g = together(f)

        if not g.is_fraction(*syms):
            return f
        else:
            p, q = g.as_numer_denom()

            p = p.expand()
            q = q.expand()

            syms = syms or None

            g = gcd(p, q, syms)

            if g is not S.One:
                p = quo(p, g, syms)
                q = quo(q, g, syms)

            return p / q

def trim(f, *syms, **flags):
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

       >>> from sympy import *

       >>> x,y = symbols('xy')
       >>> f = Function('f')

       >>> trim((f(x)**2+f(x))/f(x))
       1 + f(x)

       >>> trim((x**2+x)/x)
       1 + x

       Recursively simplify expressions:

       >>> trim(sin((f(x)**2+f(x))/f(x)))
       sin(1 + f(x))

    """
    if isinstance(f, Relational):
        return Relational(trim(f.lhs, *syms, **flags),
            trim(f.rhs, *syms, **flags), f.rel_op)
    elif isinstance(f, Matrix):
        return f.applyfunc(lambda g: trim(g, *syms, **flags))
    else:
        recursive = flags.get('recursive', True)

        def is_functional(g):
            return not (g.is_Atom or g.is_number) \
                and (not syms or g.has(*syms))

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
                        base = trim(g.base, *syms, **flags)
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
                            h = g = base**trim(g.exp, *syms, **flags)
                        else:
                            h = g = base**g.exp

                    if is_functional(h):
                        result.add(h)
                else:
                    if not recursive:
                        result.add(g)
                    else:
                        g = g.__class__(*[trim(h, *syms,
                            **flags) for h in g.args])

                        if is_functional(g):
                            result.add(g)

            return g, result

        f = sympify(f)

        if f.is_number or (syms and not f.has(*syms)):
            return f

        f = together(f.expand())
        f, terms = components(f)

        if not terms:
            return cancel(f, *syms)
        else:
            mapping, reverse = {}, {}

            for g in terms:
                mapping[g] = Temporary()
                reverse[mapping[g]] = g

            p, q = f.as_numer_denom()
            f = p.expand()/q.expand()

            H = cancel(f.subs(mapping), *syms)

            if not flags.get('extract', True):
                return H.subs(reverse)
            else:
                syms = syms or None

                def extract(f):
                    p = f.args[0]

                    for q in f.args[1:]:
                        p = gcd(p, q, syms)

                        if p.is_number:
                            return S.One, f

                    return p, Add(*[quo(g, p, syms) for g in f.args])

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

def apart(f, z, **flags):
    """Compute partial fraction decomposition of a rational function.

       Given an rational function, performing only gcd operations over
       the algebraic closue of its field of definition, compute full
       partial fraction decomposition with fractions having linear
       denominators.

       For all other kinds of expressions the input is returned in an
       unchanged form. Note however, that 'apart' function can thread
       over sums and relational operators.

       Note that no factorization of the initial denominator is needed,
       howevert at some point root finding algorithms are executed but
       for polynomials of much lower degree than the denominator.

       If given an additional flag 'formal', then even root finding
       algorithms are avoided and the result is formed as a combination
       of formal summations over implicit roots of some polynomials.

       >>> from sympy import *
       >>> x,y = symbols('xy')

       >>> apart(1/(x+2)/(x+1), x)
       1/(1 + x) - 1/(2 + x)

       >>> apart(Eq((x+1)/(x-1), E/x), x)
       1 - 2/(1 - x) == E/x

       For more information on the implemented algorithm refer to:

       [1] M. Bronstein, B. Salvy, Full partial fraction decomposition
           of rational functions, in: M. Bronstein, ed., Proceedings
           ISSAC '93, ACM Press, Kiev, Ukraine, 1993, pp. 157-160.

    """
    f = sympify(f)

    if f.is_Add:
        return Add(*[ apart(g, z, **flags) for g in f ])
    elif isinstance(f, Relational):
        return Relational(apart(f.lhs, z, **flags),
            apart(f.rhs, z, **flags), f.rel_op)
    else:
        if not f.has(z):
            return f

        if f.is_fraction(z):
            f = cancel(f, z)
        else:
            return f

        P, Q = f.as_numer_denom()

        if not Q.has(z):
            return f

        partial, r = div(P, Q, z)
        f, q, U = r / Q, Q, []

        u = Function('u')(z)
        a = Symbol('a', dummy=True)

        for k, d in enumerate(factor_.sqf(q, z)):
            n, d = k + 1, d.as_basic()
            U += [ u.diff(z, k) ]

            h = together(cancel(f * d**n, z) / u**n)

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

                g, B, _ = egcd(Q, D, z)
                b = rem(P * B / g, D, z)

                denom = (z - a)**(n-j)

                partial += RootSum(lambda r:
                    b.subs(z, r)/denom.subs(a, r), D, z)

        return partial
