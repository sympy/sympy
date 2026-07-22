from __future__ import annotations
from sympy.core import cacheit, Dummy, Ne, Integer, Rational, S, Wild
from sympy.functions import binomial, sin, cos, Piecewise, Abs
from .integrals import integrate


# TODO sin(a*x)*cos(b*x) -> sin((a+b)x) + sin((a-b)x) ?

# creating, each time, Wild's and sin/cos/Mul is expensive. Also, our match &
# subs are very slow when not cached, and if we create Wild each time, we
# effectively block caching.
#
# so we cache the pattern

# need to use a function instead of lamda since hash of lambda changes on
# each call to _pat_sincos
def _integer_instance(n):
    return isinstance(n, Integer)

@cacheit
def _pat_sincos(x):
    a = Wild('a', exclude=[x])
    n, m = [Wild(s, exclude=[x], properties=[_integer_instance])
                for s in 'nm']
    pat = sin(a*x)**n * cos(a*x)**m
    return pat, a, n, m

_u = Dummy('u')


def trigintegrate(f, x, conds='piecewise'):
    """
    Integrate f = Mul(trig) over x.

    Examples
    ========

    >>> from sympy import sin, cos, tan, sec
    >>> from sympy.integrals.trigonometry import trigintegrate
    >>> from sympy.abc import x

    >>> trigintegrate(sin(x)*cos(x), x)
    sin(x)**2/2

    >>> trigintegrate(sin(x)**2, x)
    x/2 - sin(x)*cos(x)/2

    >>> trigintegrate(tan(x)*sec(x), x)
    1/cos(x)

    >>> trigintegrate(sin(x)*tan(x), x)
    -log(sin(x) - 1)/2 + log(sin(x) + 1)/2 - sin(x)

    References
    ==========

    .. [1] https://en.wikibooks.org/wiki/Calculus/Integration_techniques

    See Also
    ========

    sympy.integrals.integrals.Integral.doit
    sympy.integrals.integrals.Integral
    """
    pat, a, n, m = _pat_sincos(x)

    f = f.rewrite('sincos')
    M = f.match(pat)

    if M is None:
        return _trigintegrate_expand(f, x)

    n, m = M[n], M[m]
    if n.is_zero and m.is_zero:
        return x
    zz = x if n.is_zero else S.Zero

    a = M[a]

    if n.is_odd or m.is_odd:
        u = _u
        n_, m_ = n.is_odd, m.is_odd

        # take smallest n or m -- to choose simplest substitution
        if n_ and m_:

            # Make sure to choose the positive one
            # otherwise an incorrect integral can occur.
            if n < 0 and m > 0:
                m_ = True
                n_ = False
            elif m < 0 and n > 0:
                n_ = True
                m_ = False
            # Both are negative so choose the smallest n or m
            # in absolute value for simplest substitution.
            elif (n < 0 and m < 0):
                n_ = n > m
                m_ = not (n > m)

            # Both n and m are odd and positive
            else:
                n_ = (n < m)      # NB: careful here, one of the
                m_ = not (n < m)  # conditions *must* be true

        #  n      m       u=C        (n-1)/2    m
        # S(x) * C(x) dx  --> -(1-u^2)       * u  du
        if n_:
            ff = -(1 - u**2)**((n - 1)/2) * u**m
            uu = cos(a*x)

        #  n      m       u=S   n         (m-1)/2
        # S(x) * C(x) dx  -->  u  * (1-u^2)       du
        elif m_:
            ff = u**n * (1 - u**2)**((m - 1)/2)
            uu = sin(a*x)

        fi = integrate(ff, u)  # XXX cyclic deps
        fx = fi.subs(u, uu)
        if conds == 'piecewise':
            return Piecewise((fx / a, Ne(a, 0)), (zz, True))
        return fx / a

    # n & m are both even
    #
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
    n_ = (Abs(n) > Abs(m))
    m_ = (Abs(m) > Abs(n))
    res = S.Zero

    if n_:
        #  2k         2 k             i             2i
        # C   = (1 - S )  = sum(i, (-) * B(k, i) * S  )
        if m > 0:
            for i in range(0, m//2 + 1):
                res += (S.NegativeOne**i * binomial(m//2, i) *
                        _sin_pow_integrate(n + 2*i, x))

        elif m == 0:
            res = _sin_pow_integrate(n, x)
        else:

            # m < 0 , |n| > |m|
            #  /
            # |
            # |    m       n
            # | cos (x) sin (x) dx =
            # |
            # |
            #/
            #                                      /
            #                                     |
            #   -1        m+1     n-1     n - 1   |     m+2     n-2
            # ________ cos (x) sin (x) + _______  |  cos (x) sin (x) dx
            #                                     |
            #   m + 1                     m + 1   |
            #                                    /

            res = (Rational(-1, m + 1) * cos(x)**(m + 1) * sin(x)**(n - 1) +
                   Rational(n - 1, m + 1) *
                   trigintegrate(cos(x)**(m + 2)*sin(x)**(n - 2), x))

    elif m_:
        #  2k         2 k            i             2i
        # S   = (1 - C ) = sum(i, (-) * B(k, i) * C  )
        if n > 0:

            #      /                            /
            #     |                            |
            #     |    m       n               |    -m         n
            #     | cos (x)*sin (x) dx  or     | cos (x) * sin (x) dx
            #     |                            |
            #    /                            /
            #
            #    |m| > |n| ; m, n >0 ; m, n belong to Z - {0}
            #       n                                         2
            #    sin (x) term is expanded here in terms of cos (x),
            #    and then integrated.
            #

            for i in range(0, n//2 + 1):
                res += (S.NegativeOne**i * binomial(n//2, i) *
                        _cos_pow_integrate(m + 2*i, x))

        elif n == 0:

            #   /
            #  |
            #  |  1
            #  | _ _ _
            #  |    m
            #  | cos (x)
            # /
            #

            res = _cos_pow_integrate(m, x)
        else:

            # n < 0 , |m| > |n|
            #  /
            # |
            # |    m       n
            # | cos (x) sin (x) dx =
            # |
            # |
            #/
            #                                      /
            #                                     |
            #    1        m-1     n+1     m - 1   |     m-2     n+2
            #  _______ cos (x) sin (x) + _______  |  cos (x) sin (x) dx
            #                                     |
            #   n + 1                     n + 1   |
            #                                    /

            res = (Rational(1, n + 1) * cos(x)**(m - 1)*sin(x)**(n + 1) +
                   Rational(m - 1, n + 1) *
                   trigintegrate(cos(x)**(m - 2)*sin(x)**(n + 2), x))

    else:
        if m == n:
            ##Substitute sin(2x)/2 for sin(x)cos(x) and then Integrate.
            res = integrate((sin(2*x)*S.Half)**m, x)
        elif (m == -n):
            if n < 0:
                # Same as the scheme described above.
                # the function argument to integrate in the end will
                # be 1, this cannot be integrated by trigintegrate.
                # Hence use sympy.integrals.integrate.
                res = (Rational(1, n + 1) * cos(x)**(m - 1) * sin(x)**(n + 1) +
                       Rational(m - 1, n + 1) *
                       integrate(cos(x)**(m - 2) * sin(x)**(n + 2), x))
            else:
                res = (Rational(-1, m + 1) * cos(x)**(m + 1) * sin(x)**(n - 1) +
                       Rational(n - 1, m + 1) *
                       integrate(cos(x)**(m + 2)*sin(x)**(n - 2), x))
    if conds == 'piecewise':
        return Piecewise((res.subs(x, a*x) / a, Ne(a, 0)), (zz, True))
    return res.subs(x, a*x) / a

def _get_linear_arg(term, x):
    """Return (a, b) such that term == a*x + b, else None."""
    a = Wild('a', exclude=[x])
    b = Wild('b', exclude=[x])
    m = term.match(a*x + b)
    if m is None:
        return None
    return m[a], m[b]

def _decompose_trig_product(expr, x):
    factors = expr.as_ordered_factors()
    result = []
    for f in factors:
        base, exp = f.as_base_exp()
        if not exp.is_Integer:
            return None
        if base.func not in (sin, cos):
            if f.is_number:
                continue
            return None
        arg = _get_linear_arg(base.args[0], x)
        if arg is None:
            return None
        a, b = arg
        if not a.is_number or a == 0:   # <-- NEW: require concrete, nonzero coefficient
            return None
        result.append((base.func, a, b, exp))
    return result

def _integrate_single_trig_term(term, x):
    """
    Directly integrate coeff*sin(a*x+b), coeff*cos(a*x+b), or a
    constant (no x), without going back through the general
    integrate() dispatcher (avoids recursion back into trigintegrate).
    Returns None if term isn't one of these simple shapes.
    """
    if not term.has(x):
        return term * x

    a = Wild('a', exclude=[x])
    b = Wild('b', exclude=[x])
    c = Wild('c', exclude=[x])

    m = term.match(c * sin(a*x + b))
    if m is not None and m[a] != 0:
        return -m[c] * cos(m[a]*x + m[b]) / m[a]

    m = term.match(c * cos(a*x + b))
    if m is not None and m[a] != 0:
        return m[c] * sin(m[a]*x + m[b]) / m[a]

    return None

def _trigintegrate_expand(expr, x):
    from sympy.simplify.fu import TR8

    decomposition = _decompose_trig_product(expr, x)
    if decomposition is None:
        return None
    if any(exp < 0 for (_, _, _, exp) in decomposition):
        return None

    total_degree = sum(abs(exp) for (_, _, _, exp) in decomposition)
    if total_degree > 12:
        return None

    expanded = TR8(expr)
    terms = expanded.as_ordered_terms()

    result = S.Zero
    for term in terms:
        integrated = _integrate_single_trig_term(term, x)
        if integrated is None:
            return None
        result += integrated
    return result

def _sin_pow_integrate(n, x):
    if n > 0:
        if n == 1:
            #Recursion break
            return -cos(x)

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

        return (Rational(-1, n) * cos(x) * sin(x)**(n - 1) +
                Rational(n - 1, n) * _sin_pow_integrate(n - 2, x))

    if n < 0:
        if n == -1:
            ##Make sure this does not come back here again.
            ##Recursion breaks here or at n==0.
            return trigintegrate(1/sin(x), x)

        # n < 0
        #  /                                                 /
        # |                                                 |
        # |    n            1               n+1     n + 2   |     n+2
        # | sin (x) dx = _______ cos (x) sin (x) + _______  |  sin (x) dx
        # |                                                 |
        # |               n + 1                     n + 1   |
        #/                                                 /
        #

        return (Rational(1, n + 1) * cos(x) * sin(x)**(n + 1) +
                Rational(n + 2, n + 1) * _sin_pow_integrate(n + 2, x))

    else:
        #n == 0
        #Recursion break.
        return x


def _cos_pow_integrate(n, x):
    if n > 0:
        if n == 1:
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

        return (Rational(1, n) * sin(x) * cos(x)**(n - 1) +
                Rational(n - 1, n) * _cos_pow_integrate(n - 2, x))

    if n < 0:
        if n == -1:
            ##Recursion break
            return trigintegrate(1/cos(x), x)

        # n < 0
        #  /                                                 /
        # |                                                 |
        # |    n            -1              n+1     n + 2   |     n+2
        # | cos (x) dx = _______ sin (x) cos (x) + _______  |  cos (x) dx
        # |                                                 |
        # |               n + 1                     n + 1   |
        #/                                                 /
        #

        return (Rational(-1, n + 1) * sin(x) * cos(x)**(n + 1) +
                Rational(n + 2, n + 1) * _cos_pow_integrate(n + 2, x))
    else:
        # n == 0
        #Recursion Break.
        return x
