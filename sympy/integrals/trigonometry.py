# -*- coding: utf-8 -*-

import sympy
from sympy.core import Symbol, Wild, S
from sympy.core.numbers import Rational
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
       >>> from sympy.abc import x

       >>> trigintegrate(sin(x)*cos(x), x)
       sin(x)**2/2

       >>> trigintegrate(sin(x)**2, x)
       x/2 - cos(x)*sin(x)/2

       http://en.wikibooks.org/wiki/Calculus/Further_integration_techniques
    """

    pat, a,n,m = _pat_sincos(x)
    ##m - cos
    ##n - sin

    M = f.match(pat)

    if M is None:
        return

    n, m = M[n], M[m]   # should always be there
    if n is S.Zero and m is S.Zero:
        return x

    a = M[a]

    if n.is_integer and m.is_integer:

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
            n_ =  (abs(n) > abs(m))
            m_ =  (abs(m) > abs(n))
            res = S.Zero

            if n_:
                #  2k       2 k             i            2i
                # C   = (1-S )  = sum(i, (-) * B(k,i) * S  )
                if m > 0 :
                    for i in range(0,m/2+1):
                        res += (-1)**i * binomial(m/2,i) * sin_pow_integrate(n+2*i, x)

                elif m == 0:
                    res=sin_pow_integrate(n,x)
                else:
                    # m < 0 , |n| > |m|
                    #  /                                                           /
                    # |                                                           |
                    # |    m       n            -1        m+1     n-1     n - 1   |     m+2     n-2
                    # | cos (x) sin (x) dx =  ________ cos (x) sin (x) + _______  |  cos (x) sin (x) dx
                    # |                                                           |
                    # |                         m + 1                     m + 1   |
                    #/                                                           /
                    #
                    #
                    res=Rational(-1,m+1)*cos(x)**(m+1)*sin(x)**(n-1) + Rational(n-1,m+1)*trigintegrate(cos(x)**(m+2)*sin(x)**(n-2),x)


            elif m_:
                #  2k        2 k            i            2i
                # S   = (1 -C ) = sum(i, (-) * B(k,i) * C  )
                if n > 0:
                    #      /                            /
                    #     |                            |
                    #     |    m       n               |    -m         n
                    #     | cos (x)*sin (x) dx  or     | cos (x) * sin (x) dx
                    #     |                            |
                    #    /                            /
                    #
                    #    |m| > |n| ; m,n >0 ; m,n belong to Z - {0}
                    #       n                                        2
                    #    sin (x) term is expanded here interms of cos (x), and then integrated.
                    for i in range(0,n/2+1):
                        res += (-1)**i * binomial(n/2,i) * cos_pow_integrate(m+2*i, x)

                elif n == 0 :
                    ##  /
                    ## |
                    #  |  1
                    #  | _ _ _
                    #  |    m
                    #  | cos (x)
                    # /
                    res= cos_pow_integrate(m,x)
                else:
                    # n < 0 , |m| > |n|
                    #  /                                                         /
                    # |                                                         |
                    # |    m       n           1        m-1     n+1     m - 1   |     m-2     n+2
                    # | cos (x) sin (x) dx = _______ cos (x) sin (x) + _______  |  cos (x) sin (x) dx
                    # |                                                         |
                    # |                       n + 1                     n + 1   |
                    #/                                                         /
                    #
                    #
                    res= Rational(1,(n+1))*cos(x)**(m-1)*sin(x)**(n+1) + Rational(m-1,n+1)*trigintegrate(cos(x)**(m-2)*sin(x)**(n+2),x)

            else :
                if m == n:
                    ##Substitute sin(2x)/2 for sin(x)cos(x) and then Integrate.
                    res=sympy.integrals.integrate((Rational(1,2)*sin(2*x))**m,x)
                elif (m == -n):
                    if n < 0:
                        ##Same as the scheme described above.
                        res= Rational(1,(n+1))*cos(x)**(m-1)*sin(x)**(n+1) + Rational(m-1,n+1)*sympy.integrals.integrate(cos(x)**(m-2)*sin(x)**(n+2),x) ##the function argument to integrate in the end will be 1 , this cannot be integrated by trigintegrate. Hence use sympy.integrals.integrate.
                    else:
                        res=Rational(-1,m+1)*cos(x)**(m+1)*sin(x)**(n-1) + Rational(n-1,m+1)*sympy.integrals.integrate(cos(x)**(m+2)*sin(x)**(n-2),x)
            return res.subs(x, a*x) / a

def sin_pow_integrate(n,x):
    if n > 0 :
        if n == 1:
            #Recursion break
            return -cos(x)
        #
        # n > 0
        #  /                                                 /
        # |                                                 |
        # |    n           -1               n-1     n - 1   |     n-2
        # | sin (x) dx =  ______ cos (x) sin (x) + _______  |  sin (x) dx
        # |                                                 |
        # |                 n                         n     |
        #/                                                 /
        #
        #

        return Rational(-1,n)*cos(x)*sin(x)**(n-1)+Rational(n-1,n)*sin_pow_integrate(n-2,x)

    if n < 0:
        if n == -1:
            ##Make sure this does not come back here again.
            ##Recursion breaks here or at n==0.
            return trigintegrate(1/sin(x),x)
        #
        # n < 0
        #  /                                                 /
        # |                                                 |
        # |    n            1               n+1     n + 2   |     n+2
        # | sin (x) dx = _______ cos (x) sin (x) + _______  |  sin (x) dx
        # |                                                 |
        # |               n + 1                     n + 1   |
        #/                                                 /
        #
        #
        return Rational(1,n+1)*cos(x)*sin(x)**(n+1) + Rational(n+2,n+1) * sin_pow_integrate(n+2,x)

    else:
        #n == 0
        #Recursion break.
        return x

def cos_pow_integrate(n,x):
    if n > 0 :
        if n==1:
            #Recursion break.
            return sin(x)

        # n > 0
        #  /                                                 /
        # |                                                 |
        # |    n            1               n-1     n - 1   |     n-2
        # | sin (x) dx =  ______ sin (x) cos (x) + _______  |  cos (x) dx
        # |                                                 |
        # |                 n                         n     |
        #/                                                 /
        #
        #

        return Rational(1,n)*sin(x)*cos(x)**(n-1)+Rational(n-1,n)*cos_pow_integrate(n-2,x)

    if n < 0:
        if n == -1:
            ##Recursion break
            return trigintegrate(1/cos(x),x)
        #
        # n < 0
        #  /                                                 /
        # |                                                 |
        # |    n            -1              n+1     n + 2   |     n+2
        # | cos (x) dx = _______ sin (x) cos (x) + _______  |  cos (x) dx
        # |                                                 |
        # |               n + 1                     n + 1   |
        #/                                                 /
        #
        #


        return Rational(-1,n+1)*sin(x)*cos(x)**(n+1) + Rational(n+2,n+1) * cos_pow_integrate(n+2,x)
    else :
        # n == 0
        #Recursion Break.
        return x



