# -*- coding: utf-8 -*-

import sympy
from sympy.core import Symbol, Wild, S
from sympy.functions import sin, cos, binomial
from sympy.core.cache import cacheit

# TODO add support for tan^m(x) * sec^n(x)
# TODO sin(a*x)*cos(b*x) -> sin((a+b)x) + sin((a-b)x) ?

# creating each time Wild's and sin/cos/Mul is expensive. Also, our match &
# subs are very slow when not cached, and if we create Wild each time, we
# effectively block caching.
#
# so we cache the pattern
@cacheit
def _pat_sincos(x):
    a, n, m = [Wild(s, exclude=[x]) for s in 'anm']
    pat = sin(a*x)**n * cos(a*x)**m

    return pat, a,n,m

_u = Symbol('u', dummy=True)



def trigintegrate(f, x):
    """Integrate f = Mul(trig) over x

       >>> from sympy import Symbol, sin, cos
       >>> from sympy.integrals.trigonometry import trigintegrate
       >>> x = Symbol('x')

       >>> trigintegrate(sin(x)*cos(x), x)
       1/2*sin(x)**2

       >>> trigintegrate(sin(x)**2, x)
       x/2 - cos(x)*sin(x)/2

       http://en.wikibooks.org/wiki/Calculus/Further_integration_techniques
    """

    pat, a,n,m = _pat_sincos(x)

    M = f.match(pat)
    if M is None:
        return

    n, m = M[n], M[m]   # should always be there
    if n is S.Zero and m is S.Zero:
        return x

    a = M[a]

    if n.is_integer and n.is_integer:

        if n.is_odd or m.is_odd:
            u = _u
            n_, m_ = n.is_odd, m.is_odd

            # take smallest n or m -- to choose simplest substitution
            if n_ and m_:
               n_ = n_ and     (n < m)  # NB: careful here, one of the
               m_ = m_ and not (n < m)  #     conditions *must* be true

            #  n      m       u=C        (n-1)/2    m
            # S(x) * C(x) dx  --> -(1-u^2)       * u  du
            if n_:
                ff = -(1-u**2)**((n-1)/2) * u**m
                uu = cos(a*x)

            #  n      m       u=S   n         (m-1)/2
            # S(x) * C(x) dx  -->  u  * (1-u^2)       du
            elif m_:
                ff = u**n * (1-u**2)**((m-1)/2)
                uu = sin(a*x)

            fi= sympy.integrals.integrate(ff, u)    # XXX cyclic deps
            fx= fi.subs(u, uu)
            return fx / a

        # n & m are even
        else:
            #               2k      2m                         2l       2l
            # we transform S (x) * C (x) into terms with only S (x) or C (x)
            #
            # example:
            #  100     4       100        2    2    100          4         2
            # S (x) * C (x) = S (x) * (1-S (x))  = S (x) * (1 + S (x) - 2*S (x))
            #
            #                  104       102     100
            #               = S (x) - 2*S (x) + S (x)
            #       2k
            # then S   is integrated with recursive formula

            # take largest n or m -- to choose simplest substitution
            n_ =     (n > m)    # NB: careful here, one of the
            m_ = not (n > m)    #     conditions *must* be true

            res = S.Zero

            if n_:
                #  2k       2 k             i            2i
                # C   = (1-S )  = sum(i, (-) * B(k,i) * S  )
                for i in range(0,m/2+1):
                    res += (-1)**i * binomial(m/2,i) * Sin_2k_integrate(n/2+i, x)

            elif m_:
                #  2k        2 k            i            2i
                # S   = (1 -C ) = sum(i, (-) * B(k,i) * C  )
                for i in range(0,n/2+1):
                    res += (-1)**i * binomial(n/2,i) * Cos_2k_integrate(m/2+i, x)

            return res.subs(x, a*x) / a

            # XXX maybe use another scheme? subst:
            #
            #  2                         2
            # S(x) = 1/2 * (1-C(2x))    C(x) = 1/2 * (1+C(2x))




# NB: PolynomialSequence in fact should be named FunctionSequence (it can
# handly all functions, not neccessary polynomials)
from sympy.functions.special.polynomials import PolynomialSequence, _x
from sympy.utilities.memoization import recurrence_memo

class Sin_2k_integrate(PolynomialSequence):
    """efficient integrate(sin(x)**2k, x)"""

    @staticmethod
    @recurrence_memo([_x])
    def calc(k, prev):
        """recursively calculate \int(sin(x)**2k, x)

        formula used:

        ⌠             n-1                ⌠
        ⎮ n          S (x)*C(x)    n-1   ⎮  n-2
        ⎮S (x)  = - ──────────── + ─── * ⎮ S (x)
        ⌡                 n         n    ⌡

        see: http://en.wikipedia.org/wiki/List_of_integrals_of_trigonometric_functions

                           n-1
        XXX maybe combine S (x)*C(x)  ->  S(n*x) + ...?
        """
        n = 2*k
        return -(sin(_x))**(n-1) * cos(_x) / n  +  prev[k-1] * (n-1)/n


class Cos_2k_integrate(PolynomialSequence):
    """efficient integrate(cos(x)**2k, x)"""

    @staticmethod
    @recurrence_memo([_x])
    def calc(k, prev):
        """recursively calculate \int(cos(x)**2k, x)

        formula used:

        ⌠             n-1                ⌠
        ⎮ n          C (x)*S(x)    n-1   ⎮  n-2
        ⎮C (x)  =   ──────────── + ─── * ⎮ C (x)
        ⌡                 n         n    ⌡

        see: http://en.wikipedia.org/wiki/List_of_integrals_of_trigonometric_functions

                           n-1
        XXX maybe combine C (x)*S(x)  ->  C(n*x) + ...?
        """
        n = 2*k
        return (cos(_x))**(n-1) * sin(_x) / n  +  prev[k-1] * (n-1)/n
