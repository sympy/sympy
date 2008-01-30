"""Module 'rewrite.py' contains advanced term rewriting methods concerning
   partial fraction decomposition, combinig together and collecting terms.
"""

from sympy.core import Basic, S, C, Symbol

from sympy.polynomials import factor_, div, quo, rem, gcd, egcd
from sympy.simplify import together

def cancel(f, *syms):
    """Cancel common factors in the given rational function.

       Given a quotient of polynomials, performing only gcd and quo
       operations in polynomial algebra, return rational function
       with numerator and denominator of minimal total degree.

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
    f = Basic.sympify(f)

    if syms and not f.has(*syms):
        return f
    elif isinstance(f, C.Add):
        return C.Add(*[ cancel(g, *syms) for g in f.args ])
    elif isinstance(f, C.Relational):
        return C.Relational(cancel(f.lhs, *syms),
            cancel(f.rhs, *syms), f.rel_op)
    else:
        g = together(f)

        if not g.is_fraction(*syms):
            return f
        else:
            p, q = g.as_numer_denom()

            syms = syms or None

            g = gcd(p, q, syms)

            if g is not S.One:
                p = quo(p, g, syms)
                q = quo(q, g, syms)

            return p / q

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

       >>> apart((x+1)/(x-1) == E/x, x)
       1 - 2/(1 - x) == E/x

       For more information on the implemented algorithm refer to:

       [1] M. Bronstein, B. Salvy, Full partial fraction decomposition
           of rational functions, in: M. Bronstein, ed., Proceedings
           ISSAC '93, ACM Press, Kiev, Ukraine, 1993, pp. 157-160.

    """
    f = Basic.sympify(f)

    if isinstance(f, C.Add):
        return C.Add(*[ apart(g, z, **flags) for g in f ])
    elif isinstance(f, C.Relational):
        return C.Relational(apart(f.lhs, z, **flags),
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

        u = C.Function('u')(z)
        A = Symbol('a', dummy=True)

        formal = flags.get('formal', False)

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

                denom = (z - A)**(n-j)

                rootof = C.RootOf(D, z)

                if not formal:
                    for root in rootof.roots():
                        partial += cancel(b.subs(z, root)) \
                            / denom.subs(A, root)
                else:
                    partial += C.Sum(b.subs(z, A) / denom, (A, rootof))

        return partial
