from calculus import defun

#----------------------------------------------------------------------------#
#                                Differentiation                             #
#----------------------------------------------------------------------------#

@defun
def difference_delta(ctx, s, n):
    r"""
    Given a sequence `(s_k)` containing at least `n+1` items, returns the
    `n`-th forward difference,

    .. math ::

        \Delta^n = \sum_{k=0}^{\infty} (-1)^{k+n} {n \choose k} s_k.
    """
    n = int(n)
    d = ctx.zero
    b = (-1) ** (n & 1)
    for k in xrange(n+1):
        d += b * s[k]
        b = (b * (k-n)) // (k+1)
    return d

@defun
def diff(ctx, f, x, n=1, method='step', scale=1, direction=0):
    r"""
    Numerically computes the derivative of `f`, `f'(x)`. Optionally,
    computes the `n`-th derivative, `f^{(n)}(x)`, for any order `n`.

    **Basic examples**

    Derivatives of a simple function::

        >>> from mpmath import *
        >>> mp.dps = 15; mp.pretty = True
        >>> diff(lambda x: x**2 + x, 1.0)
        3.0
        >>> diff(lambda x: x**2 + x, 1.0, 2)
        2.0
        >>> diff(lambda x: x**2 + x, 1.0, 3)
        0.0

    The exponential function is invariant under differentiation::

        >>> nprint([diff(exp, 3, n) for n in range(5)])
        [20.0855, 20.0855, 20.0855, 20.0855, 20.0855]

    **Method**

    One of two differentiation algorithms can be chosen with the
    ``method`` keyword argument. The two options are ``'step'``,
    and ``'quad'``. The default method is ``'step'``.

    ``'step'``:

        The derivative is computed using a finite difference
        approximation, with a small step h. This requires n+1 function
        evaluations and must be performed at (n+1) times the target
        precision. Accordingly, f must support fast evaluation at high
        precision.

    ``'quad'``:

        The derivative is computed using complex
        numerical integration. This requires a larger number of function
        evaluations, but the advantage is that not much extra precision
        is required. For high order derivatives, this method may thus
        be faster if f is very expensive to evaluate at high precision.

    With ``'quad'`` the result is likely to have a small imaginary
    component even if the derivative is actually real::

        >>> diff(sqrt, 1, method='quad')    # doctest:+ELLIPSIS
        (0.5 - 9.44...e-27j)

    **Scale**

    The scale option specifies the scale of variation of f. The step
    size in the finite difference is taken to be approximately
    eps*scale. Thus, for example if `f(x) = \cos(1000 x)`, the scale
    should be set to 1/1000 and if `f(x) = \cos(x/1000)`, the scale
    should be 1000. By default, scale = 1.

    (In practice, the default scale will work even for `\cos(1000 x)` or
    `\cos(x/1000)`. Changing this parameter is a good idea if the scale
    is something *preposterous*.)

    If numerical integration is used, the radius of integration is
    taken to be equal to scale/2. Note that f must not have any
    singularities within the circle of radius scale/2 centered around
    x. If possible, a larger scale value is preferable because it
    typically makes the integration faster and more accurate.

    **Direction**

    By default, :func:`diff` uses a central difference approximation.
    This corresponds to direction=0. Alternatively, it can compute a
    left difference (direction=-1) or right difference (direction=1).
    This is useful for computing left- or right-sided derivatives
    of nonsmooth functions:

        >>> diff(abs, 0, direction=0)
        0.0
        >>> diff(abs, 0, direction=1)
        1.0
        >>> diff(abs, 0, direction=-1)
        -1.0

    More generally, if the direction is nonzero, a right difference
    is computed where the step size is multiplied by sign(direction).
    For example, with direction=+j, the derivative from the positive
    imaginary direction will be computed.

    This option only makes sense with method='step'. If integration
    is used, it is assumed that f is analytic, implying that the
    derivative is the same in all directions.

    """
    if n == 0:
        return f(ctx.convert(x))
    orig = ctx.prec
    try:
        if method == 'step':
            ctx.prec = (orig+20) * (n+1)
            h = ctx.ldexp(scale, -orig-10)
            # Applying the finite difference formula recursively n times,
            # we get a step sum weighted by a row of binomial coefficients
            # Directed: steps x, x+h, ... x+n*h
            if direction:
                h *= ctx.sign(direction)
                steps = xrange(n+1)
                norm = h**n
            # Central: steps x-n*h, x-(n-2)*h ..., x, ..., x+(n-2)*h, x+n*h
            else:
                steps = xrange(-n, n+1, 2)
                norm = (2*h)**n
            v = ctx.difference_delta([f(x+k*h) for k in steps], n)
            v = v / norm
        elif method == 'quad':
            ctx.prec += 10
            radius = ctx.mpf(scale)/2
            def g(t):
                rei = radius*ctx.expj(t)
                z = x + rei
                return f(z) / rei**n
            d = ctx.quadts(g, [0, 2*ctx.pi])
            v = d * ctx.factorial(n) / (2*ctx.pi)
        else:
            raise ValueError("unknown method: %r" % method)
    finally:
        ctx.prec = orig
    return +v

@defun
def diffs(ctx, f, x, n=None, method='step', scale=1, direction=0):
    r"""
    Returns a generator that yields the sequence of derivatives

    .. math ::

        f(x), f'(x), f''(x), \ldots, f^{(k)}(x), \ldots

    With ``method='step'``, :func:`diffs` uses only `O(k)`
    function evaluations to generate the first `k` derivatives,
    rather than the roughly `O(k^2)` evaluations
    required if one calls :func:`diff` `k` separate times.

    With `n < \infty`, the generator stops as soon as the
    `n`-th derivative has been generated. If the exact number of
    needed derivatives is known in advance, this is further
    slightly more efficient.

    **Examples**

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> nprint(list(diffs(cos, 1, 5)))
        [0.540302, -0.841471, -0.540302, 0.841471, 0.540302, -0.841471]
        >>> for i, d in zip(range(6), diffs(cos, 1)): print i, d
        ...
        0 0.54030230586814
        1 -0.841470984807897
        2 -0.54030230586814
        3 0.841470984807897
        4 0.54030230586814
        5 -0.841470984807897

    """
    if n is None:
        n = ctx.inf
    else:
        n = int(n)

    if method != 'step':
        k = 0
        while k < n:
            yield ctx.diff(f, x, k)
            k += 1
        return

    targetprec = ctx.prec

    def getvalues(m):
        callprec = ctx.prec
        try:
            ctx.prec = workprec = (targetprec+20) * (m+1)
            h = ctx.ldexp(scale, -targetprec-10)
            if direction:
                h *= ctx.sign(direction)
                y = [f(x+h*k) for k in xrange(m+1)]
                hnorm = h
            else:
                y = [f(x+h*k) for k in xrange(-m, m+1, 2)]
                hnorm = 2*h
            return y, hnorm, workprec
        finally:
            ctx.prec = callprec

    yield f(ctx.convert(x))
    if n < 1:
        return

    if n == ctx.inf:
        A, B = 1, 2
    else:
        A, B = 1, n+1

    while 1:
        y, hnorm, workprec = getvalues(B)
        for k in xrange(A, B):
            try:
                callprec = ctx.prec
                ctx.prec = workprec
                d = ctx.difference_delta(y, k) / hnorm**k
            finally:
                ctx.prec = callprec
            yield +d
            if k >= n:
                return
        A, B = B, int(A*1.4+1)
        B = min(B, n)

@defun
def differint(ctx, f, x, n=1, x0=0):
    r"""
    Calculates the Riemann-Liouville differintegral, or fractional
    derivative, defined by

    .. math ::

        \,_{x_0}{\mathbb{D}}^n_xf(x) \frac{1}{\Gamma(m-n)} \frac{d^m}{dx^m}
        \int_{x_0}^{x}(x-t)^{m-n-1}f(t)dt

    where `f` is a given (presumably well-behaved) function,
    `x` is the evaluation point, `n` is the order, and `x_0` is
    the reference point of integration (`m` is an arbitrary
    parameter selected automatically).

    With `n = 1`, this is just the standard derivative `f'(x)`; with `n = 2`,
    the second derivative `f''(x)`, etc. With `n = -1`, it gives
    `\int_{x_0}^x f(t) dt`, with `n = -2`
    it gives `\int_{x_0}^x \left( \int_{x_0}^t f(u) du \right) dt`, etc.

    As `n` is permitted to be any number, this operator generalizes
    iterated differentiation and iterated integration to a single
    operator with a continuous order parameter.

    **Examples**

    There is an exact formula for the fractional derivative of a
    monomial `x^p`, which may be used as a reference. For example,
    the following gives a half-derivative (order 0.5)::

        >>> from mpmath import *
        >>> mp.dps = 15; mp.pretty = True
        >>> x = mpf(3); p = 2; n = 0.5
        >>> differint(lambda t: t**p, x, n)
        7.81764019044672
        >>> gamma(p+1)/gamma(p-n+1) * x**(p-n)
        7.81764019044672

    Another useful test function is the exponential function, whose
    integration / differentiation formula easy generalizes
    to arbitrary order. Here we first compute a third derivative,
    and then a triply nested integral. (The reference point `x_0`
    is set to `-\infty` to avoid nonzero endpoint terms.)::

        >>> differint(lambda x: exp(pi*x), -1.5, 3)
        0.278538406900792
        >>> exp(pi*-1.5) * pi**3
        0.278538406900792
        >>> differint(lambda x: exp(pi*x), 3.5, -3, -inf)
        1922.50563031149
        >>> exp(pi*3.5) / pi**3
        1922.50563031149

    However, for noninteger `n`, the differentiation formula for the
    exponential function must be modified to give the same result as the
    Riemann-Liouville differintegral::

        >>> x = mpf(3.5)
        >>> c = pi
        >>> n = 1+2*j
        >>> differint(lambda x: exp(c*x), x, n)
        (-123295.005390743 + 140955.117867654j)
        >>> x**(-n) * exp(c)**x * (x*c)**n * gammainc(-n, 0, x*c) / gamma(-n)
        (-123295.005390743 + 140955.117867654j)


    """
    m = max(int(ctx.ceil(ctx.re(n)))+1, 1)
    r = m-n-1
    g = lambda x: ctx.quad(lambda t: (x-t)**r * f(t), [x0, x])
    return ctx.diff(g, x, m) / ctx.gamma(m-n)

@defun
def diffun(ctx, f, n=1, **options):
    """
    Given a function f, returns a function g(x) that evaluates the nth
    derivative f^(n)(x)::

        >>> from mpmath import *
        >>> mp.dps = 15; mp.pretty = True
        >>> cos2 = diffun(sin)
        >>> sin2 = diffun(sin, 4)
        >>> cos(1.3), cos2(1.3)
        (0.267498828624587, 0.267498828624587)
        >>> sin(1.3), sin2(1.3)
        (0.963558185417193, 0.963558185417193)

    The function f must support arbitrary precision evaluation.
    See :func:`diff` for additional details and supported
    keyword options.
    """
    if n == 0:
        return f
    def g(x):
        return ctx.diff(f, x, n, **options)
    return g

@defun
def taylor(ctx, f, x, n, **options):
    r"""
    Produces a degree-`n` Taylor polynomial around the point `x` of the
    given function `f`. The coefficients are returned as a list.

        >>> from mpmath import *
        >>> mp.dps = 15; mp.pretty = True
        >>> nprint(chop(taylor(sin, 0, 5)))
        [0.0, 1.0, 0.0, -0.166667, 0.0, 0.00833333]

    The coefficients are computed using high-order numerical
    differentiation. The function must be possible to evaluate
    to arbitrary precision. See :func:`diff` for additional details
    and supported keyword options.

    Note that to evaluate the Taylor polynomial as an approximation
    of `f`, e.g. with :func:`polyval`, the coefficients must be reversed,
    and the point of the Taylor expansion must be subtracted from
    the argument:

        >>> p = taylor(exp, 2.0, 10)
        >>> polyval(p[::-1], 2.5 - 2.0)
        12.1824939606092
        >>> exp(2.5)
        12.1824939607035

    """
    return [d/ctx.factorial(i) for i, d in enumerate(ctx.diffs(f, x, n, **options))]

@defun
def pade(ctx, a, L, M):
    r"""
    Computes a Pade approximation of degree `(L, M)` to a function.
    Given at least `L+M+1` Taylor coefficients `a` approximating
    a function `A(x)`, :func:`pade` returns coefficients of
    polynomials `P, Q` satisfying

    .. math ::

        P = \sum_{k=0}^L p_k x^k

        Q = \sum_{k=0}^M q_k x^k

        Q_0 = 1

        A(x) Q(x) = P(x) + O(x^{L+M+1})

    `P(x)/Q(x)` can provide a good approximation to an analytic function
    beyond the radius of convergence of its Taylor series (example
    from G.A. Baker 'Essentials of Pade Approximants' Academic Press,
    Ch.1A)::

        >>> from mpmath import *
        >>> mp.dps = 15; mp.pretty = True
        >>> one = mpf(1)
        >>> def f(x):
        ...     return sqrt((one + 2*x)/(one + x))
        ...
        >>> a = taylor(f, 0, 6)
        >>> p, q = pade(a, 3, 3)
        >>> x = 10
        >>> polyval(p[::-1], x)/polyval(q[::-1], x)
        1.38169105566806
        >>> f(x)
        1.38169855941551

    """
    # To determine L+1 coefficients of P and M coefficients of Q
    # L+M+1 coefficients of A must be provided
    assert(len(a) >= L+M+1)

    if M == 0:
        if L == 0:
            return [ctx.one], [ctx.one]
        else:
            return a[:L+1], [ctx.one]

    # Solve first
    # a[L]*q[1] + ... + a[L-M+1]*q[M] = -a[L+1]
    # ...
    # a[L+M-1]*q[1] + ... + a[L]*q[M] = -a[L+M]
    A = ctx.matrix(M)
    for j in range(M):
        for i in range(min(M, L+j+1)):
            A[j, i] = a[L+j-i]
    v = -ctx.matrix(a[(L+1):(L+M+1)])
    x = ctx.lu_solve(A, v)
    q = [ctx.one] + list(x)
    # compute p
    p = [0]*(L+1)
    for i in range(L+1):
        s = a[i]
        for j in range(1, min(M,i) + 1):
            s += q[j]*a[i-j]
        p[i] = s
    return p, q
