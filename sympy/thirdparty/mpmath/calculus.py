"""
High-level, mostly calculus-oriented functions.

* Numerical differentiation
* Numerical polynomial operations
* Numerical root-finding
* Numerical integration

etc
"""

__docformat__ = 'plaintext'

from mptypes import *
from specfun import factorial, bernoulli_range

#----------------------------------------------------------------------------#
#                       General extrapolation methods                        #
#----------------------------------------------------------------------------#

def richardson_extrapolation(f, n, N):
    if not callable(f):
        g = f; f = lambda k: g.__getitem__(int(k))
    orig = mp.prec
    try:
        mp.prec = 2*orig
        s = mpf(0)
        for j in range(0, N+1):
            c = (n+j)**N * (-1)**(j+N) / (factorial(j) * factorial(N-j))
            s += c * f(mpf(n+j))
    finally:
        mp.prec = orig
    return +s

def shanks_extrapolation(f, n, m=1):
    if not callable(f):
        g = f; f = lambda k: g.__getitem__(int(k))
    orig = mp.prec
    try:
        mp.prec = 2*orig
        table = [f(mpf(j)) for j in range(n+m+2)]
        table2 = table[:]
        for i in range(1, m+1):
            for j in range(i, n+m+1):
                x, y, z = table[j-1], table[j], table[j+1]
                try:
                    table2[j] = (z*x - y**2) / (z + x - 2*y)
                except ZeroDivisionError:
                    return table2[j-1]
            table = table2[:]
    finally:
        mp.prec = orig
    return +table[n]

def limit(f, x, direction=-1, n=None, N=None):
    """Compute lim of f(t) as t -> x using Richardson extrapolation.
    For infinite x, the function values [f(n), ... f(n+N)] are used.
    For finite x, [f(x-direction/n)), ... f(x-direction/(n+N))] are
    used. If x is inf, f can be also be a precomputed sequence
    with a __getitem__ method."""
    if callable(f):
        if not n: n = 3 + int(mp.dps * 0.5)
        if not N: N = 2*n
    else:
        # If a sequence, take as many terms are are available
        g = f; f = lambda k: g.__getitem__(int(k))
        if not N: N = len(g)-1
        if not n: n = 0
    if   x == inf:  return richardson_extrapolation(lambda k: f(mpf(k)), n, N)
    elif x == -inf: return richardson_extrapolation(lambda k: f(mpf(-k)), n, N)
    direction *= mpf(1)
    def g(k):
        return f(x - direction/(1+k))
    return richardson_extrapolation(g, n, N)


#----------------------------------------------------------------------------#
#                                Differentiation                             #
#----------------------------------------------------------------------------#

def diff(f, x, direction=0):
    """
    Compute f'(x) using a simple finite difference approximation.

    With direction = 0, use the central difference f(x-h), f(x+h)
    With direction = 1, use the forward difference f(x), f(x+h)
    With direction = -1, use the backward difference f(x-h), f(x)

        >>> print diff(cos, 1)
        -0.841470984807897
        >>> print diff(abs, 0, 0)
        0.0
        >>> print diff(abs, 0, 1)
        1.0
        >>> print diff(abs, 0, -1)
        -1.0

    The step size is taken similar to the epsilon of the precision.
    To eliminate cancellation errors, diff temporarily doubles the
    working precision while calculating the function values.
    """
    prec = mp.prec
    extra = 5
    h = ldexp(1, -prec-extra)
    try:
        mp.prec = 2*(prec+extra)
        if   direction == 0:  return (f(x+h) - f(x-h)) * ldexp(1, prec+extra-1)
        elif direction == 1:  return (f(x+h) - f(x)) * ldexp(1, prec+extra)
        elif direction == -1: return (f(x) - f(x-h)) * ldexp(1, prec+extra)
        else:
            raise ValueError("invalid difference direction: %r" % direction)
    finally:
        mp.prec = prec

def diffc(f, x, n=1, radius=mpf(0.5)):
    """
    Compute an approximation of the nth derivative of f at the point x
    using the Cauchy integral formula. This only works for analytic
    functions. A circular path with the given radius is used.

    diffc increases the working precision slightly to avoid simple
    rounding errors. Note that, especially for large n, differentiation
    is extremely ill-conditioned, so this precaution does not
    guarantee a correct result. (Provided there are no singularities
    in the way, increasing the radius may help.)

    The returned value will be a complex number; a large imaginary part
    for a derivative that should be real may indicate a large numerical
    error.
    """
    prec = mp.prec
    try:
        mp.prec += 10
        def g(t):
            rei = radius*exp(j*t)
            z = x + rei
            return f(z) / rei**n
        d = quadts(g, 0, 2*pi)
        return d * factorial(n) / (2*pi)
    finally:
        mp.prec = prec


#----------------------------------------------------------------------------#
#                           Generic root-finding                             #
#----------------------------------------------------------------------------#

msg1 = "Cannot perform a step with the secant method because the " \
  "function values are equal at the two chosen start points. Try " \
  "different start points."

msg2 = "The derivative cannot be computed. The previous value " \
  "will be reused for the next iteration."

def secant(f, x0, x1=None, maxsteps=20, verbose=False):
    """Solve the equation f(x) = 0 using the secant method, starting
    at the given initial point x0 and performing up to `maxsteps`
    steps or quitting when the difference between successive x values
    is smaller than the epsilon of the current working precision.

    The secant method requires a second starting point x1 with both
    x0 and x1 located close to the root. If only x0 is provided, x1
    is automatically generated as x0 + 1/4."""
    weps = 2*eps
    x = x0 * mpf(1)
    if x1 is None:
        xprev = x0 + mpf(0.25)
    else:
        xprev = x1 * mpf(1)
    deriv_prev = None
    fxprev = f(xprev)
    for i in xrange(maxsteps):
        if verbose:
            print "Step", i
            print "x =", x
        fx = f(x)
        ydiff = fx - fxprev
        xdiff = x - xprev
        if verbose:
            print "f(x) =", fx
            print "xdiff = ", xdiff
            print "ydiff = ", ydiff
        try:
            deriv = xdiff / ydiff
            deriv_prev = deriv
        except ZeroDivisionError:
            if deriv_prev is None:
                raise ZeroDivisionError(msg1)
            if verbose and abs(xdiff) > weps:
                print msg2
            deriv = deriv_prev
        x, xprev = x - fx*deriv, x
        fxprev = fx
        if verbose:
            print
        if abs(xdiff) <= weps:
            break
    return x


#----------------------------------------------------------------------------#
#                                Polynomials                                 #
#----------------------------------------------------------------------------#

def polyval(coeffs, x, derivative=False):
    """
    Given coefficients [c0, c1, c2, ..., cn], evaluate
    P(x) = c0 + c1*x + c2*x**2 + ... + cn*x**n.

    If derivative=True is set, a tuple (P(x), P'(x)) is returned.
    """
    p = mpnumeric(coeffs[-1])
    q = mpf(0)
    for c in coeffs[-2::-1]:
        if derivative:
            q = p + x*q
        p = c + x*p
    if derivative:
        return p, q
    else:
        return p

def polyroots(coeffs, maxsteps=50, cleanup=True, extraprec=10, error=False):
    """
    Numerically locate all (complex) roots of a polynomial using the
    Durand-Kerner method.

    With error=True, this function returns a tuple (roots, err) where roots
    is a list of complex numbers sorted by absolute value, and err is an
    estimate of the maximum error. The polynomial should be given as a list
    of coefficients.

        >>> nprint(polyroots([24,-14,-1,1]), 4)
        [-4.0, 2.0, 3.0]
        >>> nprint(polyroots([2,3,4], error=True))
        ([(-0.375 - 0.599479j), (-0.375 + 0.599479j)], 2.22045e-16)

    """
    orig = mp.prec
    weps = +eps
    try:
        mp.prec += 10
        deg = len(coeffs) - 1
        # Must be monic
        lead = convert_lossless(coeffs[-1])
        if lead == 1:
            coeffs = map(convert_lossless, coeffs)
        else:
            coeffs = [c/lead for c in coeffs]
        f = lambda x: polyval(coeffs, x)
        roots = [mpc((0.4+0.9j)**n) for n in range(deg)]
        err = [mpf(1) for n in range(deg)]
        for step in range(maxsteps):
            if max(err).ae(0):
                break
            for i in range(deg):
                if not err[i].ae(0):
                    p = roots[i]
                    x = f(p)
                    for j in range(deg):
                        if i != j:
                            try:
                                x /= (p-roots[j])
                            except ZeroDivisionError:
                                continue
                    roots[i] = p - x
                    err[i] = abs(x)
        if cleanup:
            for i in range(deg):
                if abs(roots[i].imag) < weps:
                    roots[i] = roots[i].real
                elif abs(roots[i].real) < weps:
                    roots[i] = roots[i].imag * 1j
        roots.sort(key=lambda x: (abs(x.imag), x.real))
    finally:
        mp.prec = orig
    if error:
        err = max(err)
        err = max(err, ldexp(1, -orig+1))
        return [+r for r in roots], +err
    else:
        return [+r for r in roots]


##############################################################################
##############################################################################

#----------------------------------------------------------------------------#
#       Implementation of tanh-sinh (doubly exponential) quadrature          #
#----------------------------------------------------------------------------#

"""
The implementation of the tanh-sinh algorithm is based on the
description given in Borwein, Bailey & Girgensohn, "Experimentation
in Mathematics - Computational Paths to Discovery", A K Peters,
2003, pages 312-313.

Various documents are available online, e.g.
http://crd.lbl.gov/~dhbailey/dhbpapers/dhb-tanh-sinh.pdf
http://users.cs.dal.ca/~jborwein/tanh-sinh.pdf
"""

def transform(f, a, b):
    """Given an integrand f defined over the interval [a, b], return an
    equivalent integrand g defined on the standard interval [-1, 1].

    If a and b are finite, this is achived by means of a linear change
    of variables. If at least one point is infinite, the substitution
    t = 1/x is used."""
    a = convert_lossless(a)
    b = convert_lossless(b)
    if (a, b) == (-1, 1):
        return f
    one = mpf(1)
    half = mpf(0.5)
    # The transformation 1/x sends [1, inf] to [0, 1], which in turn
    # can be transformed to [-1, 1] the usual way. For a double
    # infinite interval, we simply evaluate the function symmetrically
    if (a, b) == (-inf, inf):
        # return transform(lambda x: (f(-1/x+1)+f(1/x-1))/x**2, 0, 1)
        def y(x):
            u = 2/(x+one)
            w = one - u
            return half * (f(w)+f(-w)) * u**2
        return y
    if a == -inf:
        # return transform(lambda x: f(-1/x+b+1)/x**2, 0, 1)
        b1 = b+1
        def y(x):
            u = 2/(x+one)
            return half * f(b1-u) * u**2
        return y
    if b == inf:
        # return transform(lambda x: f(1/x+a-1)/x**2, 0, 1)
        a1 = a-1
        def y(x):
            u = 2/(x+one)
            return half * f(a1+u) * u**2
        return y
    # Simple linear change of variables
    C = (b-a)/2
    D = (b+a)/2
    def g(x):
        return C * f(D + C*x)
    return g

def TS_estimate_error(res, prec, eps):
    """Estimate error of the calculation at the present level by
    comparing it to the results from two previous levels. The
    algorithm is that described by D. H. Bailey."""
    try:
        D1 = log(abs(res[-1]-res[-2]), 10)
        D2 = log(abs(res[-1]-res[-3]), 10)
    except ValueError:
        return eps
    D3 = -prec
    D4 = min(0, max(D1**2/D2, 2*D1, D3))
    return mpf('0.1') ** -int(D4)

def TS_guess_level(prec):
    """Guess a reasonable first level of tanh-sinh quadrature for a
    given precision. The level should not be too high, or time will
    be wasted on computing unneeded nodes, and not too low or
    the integrator may fail or have to restart unnecessarily. This
    function gives e.g.
        50 bits -> 4
        100 bits -> 5
        500 bits -> 8
        3000 bits -> 10
    These numbers are based purely on a limited amount of
    experimentation and will sometimes be wrong."""
    return int(4 + max(0, log(prec/30.0, 2)))


def TS_node(k, hn, prec, a, ar):
    """Calculate an (abscissa, weight) pair for tanh-sinh quadrature.

        x[k] = tanh(pi/2 * sinh(k*h))
        w[k] = pi/2 * cosh(k*h) / cosh(pi/2 sinh(k*h))**2
    """
    oldprec = mp.prec
    mp.prec = prec
    mp.rounding = 'up'
    t = ldexp(mpf(k), -hn)
    # We only need to calculate one exponential
    #sinht, cosht = ldexp(a-ar, -1), ldexp(a+ar, -1)
    b = exp(a - ar); br = 1/b
    sinhb, coshb = ldexp(b-br, -1), ldexp(b+br, -1)
    x, w = sinhb/coshb, (a + ar)/coshb**2
    mp.rounding = 'default'
    mp.prec = oldprec
    return x, w


TS_cache = {}

def TS_nodes(prec, m, verbose=False):
    """
    Return a list of abscissas and a list of corresponding weights for
    tanh-sinh quadrature at level m with precision prec.
    """
    if (prec, m) in TS_cache:
        return TS_cache[(prec, m)]
    eps = ldexp(1, -prec)
    t = ldexp(1, -m)
    a0 = exp(t); a0m = 1/a0
    a1 = a1m = ldexp(pi, -2)
    xs = [0]
    ws = [ldexp(pi, -1)]
    for k in xrange(1, 20 * 2**m + 1):
        a1 = a1*a0
        a1m = a1m*a0m
        x, w = TS_node(k, m, prec, a1, a1m)
        diff = abs(x-1)
        if diff <= eps:
            break
        if verbose and m > 6 and k % 300 == 150:
            # note: the number displayed is rather arbitrary. should
            # figure out how to print something that looks more like a
            # percentage
            print "Calculating nodes:", nstr(-log(diff, 10) / prec)
        xs.append(x)
        ws.append(w)
    TS_cache[(prec, m)] = (xs, ws)
    return xs, ws

def TS_eval(f, nodes, target_prec, working_prec, m, verbose=False):
    """Evaluate f at the given set of tanh-sinh nodes."""
    eps = ldexp(1, -target_prec)
    S = mpf(0)
    h = mpf(1)
    xs, ws = nodes
    res = []
    for k in xrange(1, m+1):
        if verbose:
            print "Evaluating integral (level %s of %s)" % (k, m)
        h = h / 2
        for i in xrange(0, len(xs), 2**(m-k)):
            if i % (2**(m-k+1)) != 0 or k == 1:
                if i == 0:
                    S = S + ws[0]*f(mpf(0))
                else:
                    S = S + (ws[i])*(f(-xs[i]) + f(xs[i]))
        res.append(h*S)
        if k > 2:
            err = TS_estimate_error(res, target_prec, eps)
            if verbose:
                print "Estimated error:", nstr(err)
            if err <= eps:
                break
    return +res[-1], TS_estimate_error(res, target_prec, eps)


def TS_adaptive(f, target_prec, working_prec, min_level, max_level, verbose):
    """Repeatedly attempt to integrate f, trying different levels. Quit
    as soon as the estimated error is small enough, or if that doesn't
    happen, when the max level has been tried."""
    eps = ldexp(1, -target_prec)
    for m in xrange(min_level, max_level+1):
        if verbose:
            print "Using tanh-sinh algorithm with level ", m
        nodes = TS_nodes(working_prec, m, verbose=verbose)
        s, err = TS_eval(f, nodes, target_prec, working_prec, m,
            verbose=verbose)
        steps = 2*len(nodes[0])
        if err <= eps:
            return s, err
    if verbose:
        print "Failed to reach full accuracy. Estimated error:", \
            nstr(err)
    return s, err


def quadts(f, a, b, **options):
    """
    Integrate f(x) dx over the interval [a, b], using tanh-sinh
    quadrature. Use quadts(f, (a, b), (c, d)) to calculate the
    two-dimensional integral over [a, b] x [c, d].

        >>> print quadts(lambda x: x**2, -2, 4)
        24.0
        >>> print quadts(lambda x, y: x+y, (0, 1), (0, 2))
        3.0

    Options
    =======

    target_prec
        The number of accurate bits to aim for in the result. If not
        specified, the current working precision mp.prec is used.

    working_prec
        Precision to use when evaluating the function. This should
        be set slightly higher than the target precision to eliminate
        the effects of rounding and cancellation errors.

    min_level
    max_level
        The quadts function first attempts to perform tanh-sinh
        quadrature at min_level; if that fails, at min_level+1, etc, up
        to max_level. One 'level' corresponds to roughly 2**n
        evaluation points. The levels should be integers roughly
        of the size 5-10. If not specified, reasonable values are
        inferred from the target precision.

    error
        Set to True to obtain an error estimate along with the result.

    verbose
        Set to True to display progress messages.
    """

    verbose = options.get('verbose', False)
    target_prec = options.get('target_prec', mp.prec)
    working_prec = options.get('working_prec', target_prec + 20)
    min_level = options.get('min_level', TS_guess_level(target_prec))
    max_level = options.get('max_level', min_level + 2)

    orig = mp.prec

    try:
        mp.prec = working_prec

        # Handle double integrals
        if isinstance(a, tuple):
            (a, b), (c, d) = a, b
            if c == d:
                return mpf(0)
            g = f
            # Define the inner integral to recursively call quadts. We must
            # be careful to pass along the right settings.
            def f(x):
                return quadts(lambda y: g(x,y), c, d,
                    min_level=min_level, max_level=max_level,
                    target_prec=target_prec, working_prec=working_prec)
        if a == b:
            return mpf(0)

        # Based on experience, integrals on (semi)infinite intervals require
        # a little extra work
        if inf in (abs(a), abs(b)):
            min_level += 1; max_level += 1

        # Standardize to [-1, 1] and evaluate
        f = transform(f, a, b)
        val, err = TS_adaptive(f, target_prec, working_prec,
            min_level, max_level, verbose=verbose)

        if options.get('error', False):
            return val, err
        return val
    finally:
        mp.prec = orig


def quadosc(f, a, b, period=None, zeros=None, alt=1):
    """
    Integrates f(x) from a to b where a or b is infinite and f is
    a slowly decaying oscillatory function. The zeros of f must be
    provided, either by specifying a period (suitable when f contains
    a pure sine or cosine factor) or providing a function that
    returns the nth zero (suitable when the oscillation is not
    strictly periodic).
    """
    a = convert_lossless(a)
    b = convert_lossless(b)
    if period is None and zeros is None:
        raise ValueError( \
            "either the period or zeros keyword parameter must be specified")
    if a == -inf and b == inf:
        s1 = quadosc(f, a, 0, zeros=zeros, period=period, alt=alt)
        s2 = quadosc(f, 0, b, zeros=zeros, period=period, alt=alt)
        return s1 + s2
    if a == -inf:
        if zeros:
            return quadosc(lambda x:f(-x),-b,-a, lambda n: zeros(-n), alt=alt)
        else:
            return quadosc(lambda x:f(-x),-b,-a, period=period, alt=alt)
    if b != inf:
        raise ValueError("quadosc requires an infinite integration interval")
    if not zeros:
        zeros = lambda n: n*period/2
    for n in range(1,10):
        p = zeros(n)
        if p > a:
            break
    if n >= 9:
        raise ValueError("zeros do not appear to be correctly indexed")
    if alt == 0:
        s = quadts(f, a, zeros(n+1))
        s += sumrich(lambda k: quadts(f, zeros(2*k), zeros(2*k+2)), n, inf)
    else:
        s = quadts(f, a, zeros(n))
        s += sumsh(lambda k: quadts(f, zeros(k), zeros(k+1)), n, inf)
    return s


##############################################################################
##############################################################################

#----------------------------------------------------------------------------#
#                               Numerical summation                          #
#----------------------------------------------------------------------------#

def sumrich(f, a, b, n=None, N=None):
    """Sum f(k) for k = a..b using Richardson extrapolation. This
    function is essentially equivalent to using limit() on the
    sequence of partial sums."""
    assert b == inf
    if not n: n = 3 + int(mp.dps * 0.5)
    if not N: N = 2*n
    orig = mp.prec
    try:
        mp.prec = 2*orig
        s = mpf(0)
        tbl = []
        for k in range(a, a+n+N+1):
            s += f(mpf(k))
            tbl.append(s)
        s = richardson_extrapolation(lambda k: tbl[int(k)], n, N)
    finally:
        mp.prec = orig
    return +s

def sumsh(f, a, b, n=None, m=None):
    """Sum f(k) for k = a..b using an n-term Shanks
    transformation. With m > 1, the Shanks transformation is applied
    recursively m times.

    Shanks summation often works well for slowly convergent and/or
    alternating Taylor series."""
    assert b == inf
    if not n: n = 5 + int(mp.dps * 1.2)
    if not m: m = 2 + n//3
    orig = mp.prec
    try:
        mp.prec = 2*orig
        s = mpf(0)
        tbl = []
        for k in range(a, a+n+m+2):
            s += f(mpf(k))
            tbl.append(s)
        s = shanks_extrapolation(tbl, n, m)
    finally:
        mp.prec = orig
    return +s

@extraprec(15, normalize_output=True)
def sumem(f, a=0, b=inf, N=None, integral=None, fderiv=None, error=False,
    verbose=False):
    """
    Calculate the sum of f(n) for n = a..b using Euler-Maclaurin
    summation. This algorithm is efficient for slowly convergent
    nonoscillatory sums; the essential condition is that f must be
    analytic. The method relies on approximating the sum by an
    integral, so f must be smooth and well-behaved enough to be
    integrated numerically.

    A tuple (s, err) is returned where s is the calculated sum and err
    is the estimated magnitude of the error. With verbose=True,
    detailed information about progress and errors is printed.

        >>> mp.dps = 15
        >>> s, err = sumem(lambda n: 1/n**2, 1, inf, error=True)
        >>> print s
        1.64493406684823
        >>> print pi**2 / 6
        1.64493406684823
        >>> nprint(err)
        2.22045e-16

    N is the number of terms to compute directly before using the
    Euler-Maclaurin formula to approximate the tail. It must be set
    high enough; often roughly N ~ dps is the right size.

    High-order derivatives of f are also needed. By default, these
    are computed using numerical integration, which is the most
    expensive part of the calculation. The default method assumes
    that all poles of f are located close to the origin. A custom
    nth derivative function fderiv(x, n) can be provided as a
    keyword parameter.

    This is much more efficient:

        >>> f = lambda n: 1/n**2
        >>> fp = lambda x, n: (-1)**n * factorial(n+1) * x**(-2-n)
        >>> mp.dps = 50
        >>> print sumem(lambda n: 1/n**2, 1, inf, fderiv=fp)
        1.6449340668482264364724151666460251892189499012068
        >>> print pi**2 / 6
        1.6449340668482264364724151666460251892189499012068

    If b = inf, f and its derivatives are all assumed to vanish
    at infinity. It is assumed that a is finite, so doubly
    infinite sums cannot be evaluated directly.
    """
    if N is None:
        N = 3*mp.dps + 20
    a, b, N = mpf(a), mpf(b), mpf(N)
    infinite = (b == inf)
    weps = eps * 2**8
    if verbose:
        print "Summing f(k) from k = %i to %i" % (a, a+N-1)
    S = sum(f(mpf(k)) for k in xrange(a, a+N))
    if integral is None:
        if verbose:
            print "Integrating f(x) from x = %i to %s" % (a+N, nstr(b))
        I, ierr = quadts(f, a+N, b, error=1)
    else:
        I, ierr = integral(a+N, b), mpf(0)
    # There is little hope if the tail cannot be integrated
    # accurately. Estimate magnitude of tail as the error.
    if ierr > weps:
        if verbose:
            print "Failed to converge to target accuracy (integration failed)"
        if error:
            return S+I, abs(I) + ierr
        else:
            return S+I
    if infinite:
        C = f(a+N) / 2
    else:
        C = (f(a+N) + f(b)) / 2
    # Default (inefficient) approach for derivatives
    if not fderiv:
        fderiv = lambda x, n: diffc(f, x, n, radius=N*0.75)
    k = 1
    prev = 0
    if verbose:
        print "Summing tail"
    B = bernoulli_range()
    fac = 2
    while 1:
        if infinite:
            D = fderiv(a+N, 2*k-1)
        else:
            D = fderiv(a+N, 2*k-1) - fderiv(b, 2*k-1)
        # B(2*k) / fac(2*k)
        term = B.next() / fac * D
        mag = abs(term)
        if verbose:
            print "term", k, "magnitude =", nstr(mag)
        # Error can be estimated as the magnitude of the smallest term
        if k >= 2:
            if mag < weps:
                if verbose:
                    print "Converged to target accuracy"
                res, err = I + C + S, eps * 2**15
                break
            if mag > abs(prev):
                if verbose:
                    print "Failed to converge to target accuracy (N too low)"
                res, err = I + C + S, abs(term)
                break
        S -= term
        k += 1
        fac *= (2*k) * (2*k-1)
        prev = term
    if isinstance(res, mpc) and not isinstance(I, mpc):
        res, err = res.real, err
    if error:
        return res, err
    else:
        return res

#----------------------------------------------------------------------------#
#                                  ODE solvers                               #
#----------------------------------------------------------------------------#

def smul(a, x):
    """Multiplies the vector "x" by the scalar "a"."""
    R = []
    for i in range(len(x)):
        R.append(a*x[i])
    return R

def vadd(*args):
    """Adds vectors "x", "y", ... together."""
    assert len(args) >= 2
    n = len(args[0])
    for x in args:
        assert len(x) == n
    R = []
    for i in range(n):
        s = 0.
        for x in args:
            s += x[i]
        R.append(s)
    return R

def ODE_step_euler(x, y, h, derivs):
    """
    Advances the solution y(x) from x to x+h using the Euler method.

    derivs .... a python function f(x, (y1, y2, y3, ...)) returning
    a tuple (y1', y2', y3', ...) where y1' is the derivative of y1 at x.
    """
    X = derivs(y,x)
    return vadd(y, smul(h, X))

half = mpf(0.5)

def ODE_step_rk4(x, y, h, derivs):
    """
    Advances the solution y(x) from x to x+h using the 4th-order Runge-Kutta
    method.

    derivs .... a python function f(x, (y1, y2, y3, ...)) returning
    a tuple (y1', y2', y3', ...) where y1' is the derivative of y1 at x.
    """
    h2 = h/2
    third = mpf(1)/3
    sixth = mpf(1)/6
    k1 = smul(h, derivs(y, x))
    k2 = smul(h, derivs(vadd(y, smul(half, k1)), x+h2))
    k3 = smul(h, derivs(vadd(y, smul(half, k2)), x+h2))
    k4 = smul(h, derivs(vadd(y, k3), x+h))
    return vadd(y, smul(sixth, k1), smul(third, k2), smul(third, k3), 
            smul(sixth, k4))

def odeint(derivs, x0, t_list, step=ODE_step_rk4):
    """
    Given the list t_list of values, returns the solution at these points.
    """
    x = x0
    result = [x]
    for i in range(len(t_list)-1):
        dt = t_list[i+1] - t_list[i]
        x = step(t_list[i], x, dt, derivs)
        result.append(x)
    return result

#----------------------------------------------------------------------------#
#                              Approximation methods                         #
#----------------------------------------------------------------------------#

# The Chebyshev approximation formula is given at:
# http://mathworld.wolfram.com/ChebyshevApproximationFormula.html

# The only major changes in the following code is that we return the
# expanded polynomial coefficients instead of Chebyshev coefficients,
# and that we automatically transform [a,b] -> [-1,1] and back
# for convenience.

# Coefficient in Chebyshev approximation
def chebcoeff(f,a,b,j,N):
    s = mpf(0)
    h = mpf(0.5)
    for k in range(1, N+1):
        t = cos(pi*(k-h)/N)
        s += f(t*(b-a)*h + (b+a)*h) * cos(pi*j*(k-h)/N)
    return 2*s/N

# Generate Chebyshev polynomials T_n(ax+b) in expanded form
def chebT(a=1, b=0):
    Tb = [1]
    yield Tb
    Ta = [b, a]
    while 1:
        yield Ta
        # Recurrence: T[n+1](ax+b) = 2*(ax+b)*T[n](ax+b) - T[n-1](ax+b)
        Tmp = [0] + [2*a*t for t in Ta]
        for i, c in enumerate(Ta): Tmp[i] += 2*b*c
        for i, c in enumerate(Tb): Tmp[i] -= c
        Ta, Tb = Tmp, Ta

def chebyfit(f,a,b,N,error=False):
    """Chebyshev approximation: returns coefficients of a degree N-1
    polynomial that approximates f on the interval [a, b]. With error=True,
    also returns an estimate of the maximum error."""
    orig = mp.prec
    try:
        mp.prec = orig + int(N**0.5) + 20
        c = [chebcoeff(f,a,b,k,N) for k in range(N)]
        d = [mpf(0)] * N
        d[0] = -c[0]/2
        h = mpf(0.5)
        T = chebT(mpf(2)/(b-a), mpf(-1)*(b+a)/(b-a))
        for k in range(N):
            Tk = T.next()
            for i in range(len(Tk)):
                d[i] += c[k]*Tk[i]
        # Estimate maximum error
        err = mpf(0)
        for k in range(N):
            x = cos(pi*k/N) * (b-a)*h + (b+a)*h
            err = max(err, abs(f(x) - polyval(d, x)))
    finally:
        mp.prec = orig
        if error:
            return d, +err
        else:
            return d


#----------------------------------------------------------------------------#
#                 Lattice reduction and constant recognition                 #
#----------------------------------------------------------------------------#


"""
This is a fairly direct translation to Python of the pseudocode given by
David Bailey, "The PSLQ Integer Relation Algorithm":
http://www.cecm.sfu.ca/organics/papers/bailey/paper/html/node3.html

The stopping criteria are NOT yet properly implemented.
"""
def pslq(x, eps=None):
    """
    Given a vector of real numbers x = [x1, x2, ..., xn], pslq(x) uses the
    PSLQ algorithm to find a list of integers [c1, c2, ..., cn] such that
    c1*x1 + c2*x2 + ... + cn*xn = 0 approximately.
    """
    n = len(x)
    assert n >= 1
    prec = mp.prec
    assert prec >= 53
    target = prec // max(2,n)
    if target < 30:
        if target < 5:
            print "Warning: precision for PSLQ may be too low"
        target = int(prec * 0.75)
    if eps is None:
        eps = mpf(2)**(-target)
    x = [None] + x
    g = sqrt(mpf(4)/3)
    A = {}
    B = {}
    H = {}
    # Initialization
    # step 1
    for i in range(1, n+1):
        for j in range(1, n+1):
            A[i,j] = B[i,j] = mpf(int(i == j))
            H[i,j] = mpf(0)
    # step 2
    s = [None] + [mpf(0)] * n
    for k in range(1, n+1):
        t = mpf(0)
        for j in range(k, n+1):
            t += x[j]**2
        s[k] = sqrt(t)
    t = s[1]
    y = x[:]
    for k in range(1, n+1):
        y[k] = x[k] / t
        s[k] = s[k] / t
    # step 3
    for i in range(1, n+1):
        for j in range(i+1, n): H[i,j] = mpf(0)
        if i <= n-1: H[i,i] = s[i+1]/s[i]
        for j in range(1, i): H[i,j] = -y[i]*y[j]/(s[j]*s[j+1])
    # step 4
    for i in range(2, n+1):
        for j in range(i-1, 0, -1):
            t = floor(H[i,j]/H[j,j] + 0.5)
            y[j] = y[j] + t*y[i]
            for k in range(1, j+1):
                H[i,k] = H[i,k] - t*H[j,k]
            for k in range(1, n+1):
                A[i,k] = A[i,k] - t*A[j,k]
                B[k,j] = B[k,j] + t*B[k,i]
    # Main algorithm
    for REP in range(100):
        # step 1
        m = -1
        szmax = -1
        for i in range(1, n):
            h = H[i,i]
            sz = sqrt(mpf(4)/3)**i * abs(h)
            if sz > szmax:
                m = i
                szmax = sz
        # step 2
        y[m], y[m+1] = y[m+1], y[m]
        tmp = {}
        for i in range(1,n+1): H[m,i], H[m+1,i] = H[m+1,i], H[m,i]
        for i in range(1,n+1): A[m,i], A[m+1,i] = A[m+1,i], A[m,i]
        for i in range(1,n+1): B[i,m], B[i,m+1] = B[i,m+1], B[i,m]
        # step 3
        if m <= n - 2:
            t0 = sqrt(H[m,m]**2 + H[m,m+1]**2)
            t1 = H[m,m] / t0
            t2 = H[m,m+1] / t0
            for i in range(m, n+1):
                t3 = H[i,m]
                t4 = H[i,m+1]
                H[i,m] = t1*t3+t2*t4
                H[i,m+1] = -t2*t3+t1*t4
        # step 4
        for i in range(m+1, n+1):
            for j in range(min(i-1, m+1), 0, -1):
                try:
                    t = floor(H[i,j]/H[j,j] + 0.5)
                # XXX
                except ZeroDivisionError:
                    break
                y[j] = y[j] + t*y[i]
                for k in range(1, j+1):
                    H[i,k] = H[i,k] - t*H[j,k]
                for k in range(1, n+1):
                    A[i,k] = A[i,k] - t*A[j,k]
                    B[k,j] = B[k,j] + t*B[k,i]
        for i in range(1, n+1):
            if abs(y[i]) < eps:
                vec = [int(int(B[j,i])) for j in range(1,n+1)]
                if max(abs(v) for v in vec) < 10**6:
                    return vec
    return None

def findpoly(x, n=1):
    if x == 0:
        return [0, 1]
    xs = [mpf(1)]
    for i in range(1,n+1):
        xs.append(x**i)
        a = pslq(xs)
        if a is not None:
            return a

def fracgcd(p, q):
    x, y = p, q
    while y:
        x, y = y, x % y
    if x != 1:
        p //= x
        q //= x
    if q == 1:
        return p
    return p, q

def pslqstring(r, constants):
    q = r[0]
    r = r[1:]
    s = []
    for i in range(len(r)):
        p = r[i]
        if p:
            z = fracgcd(-p,q)
            cs = constants[i][1]
            if cs == '1':
                cs = ''
            else:
                cs = '*' + cs
            if isinstance(z, (int, long)):
                if z > 0: term = str(z) + cs
                else:     term = ("(%s)" % z) + cs
            else:
                term = ("(%s/%s)" % z) + cs
            s.append(term)
    s = ' + '.join(s)
    if '+' in s or '*' in s:
        s = '(' + s + ')'
    return s or '0'

def prodstring(r, constants):
    q = r[0]
    r = r[1:]
    num = []
    den = []
    for i in range(len(r)):
        p = r[i]
        if p:
            z = fracgcd(-p,q)
            cs = constants[i][1]
            if isinstance(z, (int, long)):
                if abs(z) == 1: t = cs
                else:           t = '%s**%s' % (cs, abs(z))
                ([num,den][z<0]).append(t)
            else:
                t = '%s**(%s/%s)' % (cs, abs(z[0]), z[1])
                ([num,den][z[0]<0]).append(t)
    num = '*'.join(num)
    den = '*'.join(den)
    if num and den: return "(%s)/(%s)" % (num, den)
    if num: return num
    if den: return "1/(%s)" % den

def quadraticstring(t,a,b,c):
    if c < 0:
        a,b,c = -a,-b,-c
    u1 = (-b+sqrt(b**2-4*a*c))/(2*c)
    u2 = (-b-sqrt(b**2-4*a*c))/(2*c)
    if abs(u1-t) < abs(u2-t):
        if b:  s = '((%s+sqrt(%s))/%s)' % (-b,b**2-4*a*c,2*c)
        else:  s = '(sqrt(%s)/%s)' % (-4*a*c,2*c)
    else:
        if b:  s = '((%s-sqrt(%s))/%s)' % (-b,b**2-4*a*c,2*c)
        else:  s = '(-sqrt(%s)/%s)' % (-4*a*c,2*c)
    return s

# Transformation y = f(x,c), with inverse function x = f(y,c)
# The third entry indicates whether the transformation is
# redundant when c = 1
transforms = [
  (lambda x,c: x*c, '$y/$c', 0),
  (lambda x,c: x/c, '$c*$y', 1),
  (lambda x,c: c/x, '$c/$y', 0),
  (lambda x,c: (x*c)**2, 'sqrt($y)/$c', 0),
  (lambda x,c: (x/c)**2, '$c*sqrt($y)', 1),
  (lambda x,c: (c/x)**2, '$c/sqrt($y)', 0),
  (lambda x,c: c*x**2, 'sqrt($y)/sqrt($c)', 1),
  (lambda x,c: x**2/c, 'sqrt($c)*sqrt($y)', 1),
  (lambda x,c: c/x**2, 'sqrt($c)/sqrt($y)', 1),
  (lambda x,c: sqrt(x*c), '$y**2/$c', 0),
  (lambda x,c: sqrt(x/c), '$c*$y**2', 1),
  (lambda x,c: sqrt(c/x), '$c/$y**2', 0),
  (lambda x,c: c*sqrt(x), '$y**2/$c**2', 1),
  (lambda x,c: sqrt(x)/c, '$c**2*$y**2', 1),
  (lambda x,c: c/sqrt(x), '$c**2/$y**2', 1),
  (lambda x,c: exp(x*c), 'log($y)/$c', 0),
  (lambda x,c: exp(x/c), '$c*log($y)', 1),
  (lambda x,c: exp(c/x), '$c/log($y)', 0),
  (lambda x,c: c*exp(x), 'log($y/$c)', 1),
  (lambda x,c: exp(x)/c, 'log($c*$y)', 1),
  (lambda x,c: c/exp(x), 'log($c/$y)', 0),
  (lambda x,c: log(x*c), 'exp($y)/$c', 0),
  (lambda x,c: log(x/c), '$c*exp($y)', 1),
  (lambda x,c: log(c/x), '$c/exp($y)', 0),
  (lambda x,c: c*log(x), 'exp($y/$c)', 1),
  (lambda x,c: log(x)/c, 'exp($c*$y)', 1),
  (lambda x,c: c/log(x), 'exp($c/$y)', 0),
]

def identify(x, constants=[], full=False, maxcoeff=1000, tolerance=None,
    quadratics=True, verbose=False):
    """"
    This function attempts to find a symbolic expression for the given
    quantity x. It can identify simple algebraic numbers, as well as
    simple combinations of the given list of base constants, and
    exponentials or logarithms thereof.

    The base constants should be given as strings representing mpmath
    expressions.

    If a match is found, a mathematical formula is returned as a string.
    With full=True, a list of matching formulas is returned.

    In order not to produce spurious approximations, high precision
    should be used; preferrably 50 digits or more.

    Examples:

        >>> mp.dps = 15
        >>> identify(0.22222222222222222)
        '(2/9)'

        >>> mp.dps = 50
        >>> identify(3*pi + 4*sqrt(2), ['pi','sqrt(2)'])
        ...
        (3*pi + 4*sqrt(2))

    Further example identifications that should work (many redundant
    results may be found if run with full=True):

        mp.dps = 50
        base = ['sqrt(2)','pi','log(2)']
        identify(0.25, base)
        identify(3*pi + 2*sqrt(2) + 5*log(2)/7, base)
        identify(exp(pi+2), base)
        identify(1/(3+sqrt(2)), base)
        identify(sqrt(2)/(3*pi+4), base)
        identify(5**(mpf(1)/3)*pi*log(2)**2, base)
    """

    solutions = []

    def addsolution(s):
        if verbose: print "Found: ", s
        solutions.append(s)

    x = mpf(x)

    # Further along, x will be assumed positive
    if x == 0:
        if full: return ['0']
        else:    return '0'

    if x < 0:
        sol = identify(-x, constants, full, maxcoeff, tolerance, quadratics, verbose)
        if sol is None:
            return sol
        if full:
            return ["-(%s)"%s for s in sol]
        else:
            return "-(%s)" % sol

    sols = []
    if tolerance:
        weps = mpf(tolerance)
    else:
        weps = eps**0.7

    if isinstance(constants, dict):
        constants = [(mpf(v), name) for (name, v) in constants.items()]
    else:
        import mpmath
        constants = [(eval(p, mpmath.__dict__), p) for p in constants]

    # We always want to find at least rational terms
    if 1 not in [value for (name, value) in constants]:
        constants = [(mpf(1), '1')] + constants

    # PSLQ with simple algebraic and functional transformations
    for ft, ftn, red in transforms:
        for c, cn in constants:
            if red and cn == '1':
                continue
            t = ft(x,c)
            # Prevent exponential transforms from wreaking havoc
            if abs(t) > maxcoeff**2 or abs(t) < weps:
                continue
            # Linear combination of base constants
            r = pslq([t] + [a[0] for a in constants], weps)
            s = None
            if r is not None and max(abs(uw) for uw in r) <= maxcoeff and r[0]:
                s = pslqstring(r, constants)
            # Quadratic algebraic numbers
            elif quadratics:
                q = pslq([mpf(1), t, t**2], weps)
                if q is not None and len(q) == 3 and q[2]:
                    aa,bb,cc = q
                    if max(abs(aa),abs(bb),abs(cc)) <= maxcoeff:
                        s = quadraticstring(t,aa,bb,cc)
            if s:
                if cn == '1' and ('/$c' in ftn):
                    s = ftn.replace('$y', s).replace('/$c', '')
                else:
                    s = ftn.replace('$y', s).replace('$c', cn)
                addsolution(s)
                if not full: return solutions[0]

            if verbose:
                print "."

    # Check for a direct multiplicative formula
    if x != 1:
        # Allow fractional powers of fractions
        ilogs = [2,3,5,7]
        # Watch out for existing fractional powers of fractions
        logs = []
        for a, s in constants:
            if not sum(bool(findpoly(log(a)/log(i),1)) for i in ilogs):
                logs.append((log(a), s))
        logs = [(log(i),str(i)) for i in ilogs] + logs
        r = pslq([log(x)] + [a[0] for a in logs], weps)
        if r is not None and max(abs(uw) for uw in r) <= maxcoeff and r[0]:
            addsolution(prodstring(r, logs))
            if not full: return solutions[0]

    if full:
        return sorted(solutions, key=len)
    else:
        return None
