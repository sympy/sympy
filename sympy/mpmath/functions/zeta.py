from functions import defun, defun_wrapped, defun_static

@defun
def stieltjes(ctx, n, a=1):
    n = ctx.convert(n)
    a = ctx.convert(a)
    if n < 0:
        return ctx.bad_domain("Stieltjes constants defined for n >= 0")
    if hasattr(ctx, "stieltjes_cache"):
        stieltjes_cache = ctx.stieltjes_cache
    else:
        stieltjes_cache = ctx.stieltjes_cache = {}
    if a == 1:
        if n == 0:
            return +ctx.euler
        if n in stieltjes_cache:
            prec, s = stieltjes_cache[n]
            if prec >= ctx.prec:
                return +s
    mag = 1
    def f(x):
        xa = x/a
        v = (xa-ctx.j)*ctx.ln(a-ctx.j*x)**n/(1+xa**2)/(ctx.exp(2*ctx.pi*x)-1)
        return ctx._re(v) / mag
    orig = ctx.prec
    try:
        # Normalize integrand by approx. magnitude to
        # speed up quadrature (which uses absolute error)
        if n > 50:
            ctx.prec = 20
            mag = ctx.quad(f, [0,ctx.inf], maxdegree=3)
        ctx.prec = orig + 10 + int(n**0.5)
        s = ctx.quad(f, [0,ctx.inf], maxdegree=20)
        v = ctx.ln(a)**n/(2*a) - ctx.ln(a)**(n+1)/(n+1) + 2*s/a*mag
    finally:
        ctx.prec = orig
    if a == 1 and ctx.isint(n):
        stieltjes_cache[n] = (ctx.prec, v)
    return +v

@defun_wrapped
def siegeltheta(ctx, t):
    if ctx._im(t):
        # XXX: cancellation occurs
        a = ctx.loggamma(0.25+0.5j*t)
        b = ctx.loggamma(0.25-0.5j*t)
        return -ctx.ln(ctx.pi)/2*t - 0.5j*(a-b)
    else:
        if ctx.isinf(t):
            return t
        return ctx._im(ctx.loggamma(0.25+0.5j*t)) - ctx.ln(ctx.pi)/2*t

@defun_wrapped
def grampoint(ctx, n):
    # asymptotic expansion, from
    # http://mathworld.wolfram.com/GramPoint.html
    g = 2*ctx.pi*ctx.exp(1+ctx.lambertw((8*n+1)/(8*ctx.e)))
    return ctx.findroot(lambda t: ctx.siegeltheta(t)-ctx.pi*n, g)

@defun_wrapped
def siegelz(ctx, t):
    v = ctx.expj(ctx.siegeltheta(t))*ctx.zeta(0.5+ctx.j*t)
    if ctx._is_real_type(t):
        return ctx._re(v)
    return v

_zeta_zeros = [
14.134725142,21.022039639,25.010857580,30.424876126,32.935061588,
37.586178159,40.918719012,43.327073281,48.005150881,49.773832478,
52.970321478,56.446247697,59.347044003,60.831778525,65.112544048,
67.079810529,69.546401711,72.067157674,75.704690699,77.144840069,
79.337375020,82.910380854,84.735492981,87.425274613,88.809111208,
92.491899271,94.651344041,95.870634228,98.831194218,101.317851006,
103.725538040,105.446623052,107.168611184,111.029535543,111.874659177,
114.320220915,116.226680321,118.790782866,121.370125002,122.946829294,
124.256818554,127.516683880,129.578704200,131.087688531,133.497737203,
134.756509753,138.116042055,139.736208952,141.123707404,143.111845808,
146.000982487,147.422765343,150.053520421,150.925257612,153.024693811,
156.112909294,157.597591818,158.849988171,161.188964138,163.030709687,
165.537069188,167.184439978,169.094515416,169.911976479,173.411536520,
174.754191523,176.441434298,178.377407776,179.916484020,182.207078484,
184.874467848,185.598783678,187.228922584,189.416158656,192.026656361,
193.079726604,195.265396680,196.876481841,198.015309676,201.264751944,
202.493594514,204.189671803,205.394697202,207.906258888,209.576509717,
211.690862595,213.347919360,214.547044783,216.169538508,219.067596349,
220.714918839,221.430705555,224.007000255,224.983324670,227.421444280,
229.337413306,231.250188700,231.987235253,233.693404179,236.524229666,
]

def _load_zeta_zeros(url):
    import urllib
    d = urllib.urlopen(url)
    L = [float(x) for x in d.readlines()]
    # Sanity check
    assert round(L[0]) == 14
    _zeta_zeros[:] = L

@defun
def zetazero(ctx, n, url='http://www.dtc.umn.edu/~odlyzko/zeta_tables/zeros1'):
    n = int(n)
    if n < 0:
        return ctx.zetazero(-n).conjugate()
    if n == 0:
        raise ValueError("n must be nonzero")
    if n > len(_zeta_zeros) and n <= 100000:
        _load_zeta_zeros(url)
    if n > len(_zeta_zeros):
        raise NotImplementedError("n too large for zetazeros")
    return ctx.mpc(0.5, ctx.findroot(ctx.siegelz, _zeta_zeros[n-1]))

@defun_wrapped
def riemannr(ctx, x):
    if x == 0:
        return ctx.zero
    # Check if a simple asymptotic estimate is accurate enough
    if abs(x) > 1000:
        a = ctx.li(x)
        b = 0.5*ctx.li(ctx.sqrt(x))
        if abs(b) < abs(a)*ctx.eps:
            return a
    if abs(x) < 0.01:
        # XXX
        ctx.prec += int(-ctx.log(abs(x),2))
    # Sum Gram's series
    s = t = ctx.one
    u = ctx.ln(x)
    k = 1
    while abs(t) > abs(s)*ctx.eps:
        t = t * u / k
        s += t / (k * ctx._zeta_int(k+1))
        k += 1
    return s

@defun_static
def primepi(ctx, x):
    x = int(x)
    if x < 2:
        return 0
    return len(ctx.list_primes(x))

@defun_wrapped
def primepi2(ctx, x):
    x = int(x)
    if x < 2:
        return ctx.mpi(0,0)
    if x < 2657:
        return ctx.mpi(ctx.primepi(x))
    mid = ctx.li(x)
    # Schoenfeld's estimate for x >= 2657, assuming RH
    err = ctx.sqrt(x,rounding='u')*ctx.ln(x,rounding='u')/8/ctx.pi(rounding='d')
    a = ctx.floor((ctx.mpi(mid)-err).a, rounding='d')
    b = ctx.ceil((ctx.mpi(mid)+err).b, rounding='u')
    return ctx.mpi(a, b)

@defun_wrapped
def primezeta(ctx, s):
    if ctx.isnan(s):
        return s
    if ctx.re(s) <= 0:
        raise ValueError("prime zeta function defined only for re(s) > 0")
    if s == 1:
        return ctx.inf
    if s == 0.5:
        return ctx.mpc(ctx.ninf, ctx.pi)
    r = ctx.re(s)
    if r > ctx.prec:
        return 0.5**s
    else:
        wp = ctx.prec + int(r)
        def terms():
            orig = ctx.prec
            # zeta ~ 1+eps; need to set precision
            # to get logarithm accurately
            k = 0
            while 1:
                k += 1
                u = ctx.moebius(k)
                if not u:
                    continue
                ctx.prec = wp
                t = u*ctx.ln(ctx.zeta(k*s))/k
                if not t:
                    return
                #print ctx.prec, ctx.nstr(t)
                ctx.prec = orig
                yield t
    return ctx.sum_accurately(terms)

# TODO: for bernpoly and eulerpoly, ensure that all exact zeros are covered

@defun_wrapped
def bernpoly(ctx, n, z):
    # Slow implementation:
    #return sum(ctx.binomial(n,k)*ctx.bernoulli(k)*z**(n-k) for k in xrange(0,n+1))
    n = int(n)
    if n < 0:
        raise ValueError("Bernoulli polynomials only defined for n >= 0")
    if z == 0 or (z == 1 and n > 1):
        return ctx.bernoulli(n)
    if z == 0.5:
        return (ctx.ldexp(1,1-n)-1)*ctx.bernoulli(n)
    if n <= 3:
        if n == 0: return z ** 0
        if n == 1: return z - 0.5
        if n == 2: return (6*z*(z-1)+1)/6
        if n == 3: return z*(z*(z-1.5)+0.5)
    if abs(z) == ctx.inf:
        return z ** n
    if z != z:
        return z
    if abs(z) > 2:
        def terms():
            t = ctx.one
            yield t
            r = ctx.one/z
            k = 1
            while k <= n:
                t = t*(n+1-k)/k*r
                if not (k > 2 and k & 1):
                    yield t*ctx.bernoulli(k)
                k += 1
        return ctx.sum_accurately(terms) * z**n
    else:
        def terms():
            yield ctx.bernoulli(n)
            t = ctx.one
            k = 1
            while k <= n:
                t = t*(n+1-k)/k * z
                m = n-k
                if not (m > 2 and m & 1):
                    yield t*ctx.bernoulli(m)
                k += 1
        return ctx.sum_accurately(terms)

@defun_wrapped
def eulerpoly(ctx, n, z):
    n = int(n)
    if n < 0:
        raise ValueError("Euler polynomials only defined for n >= 0")
    if n <= 2:
        if n == 0: return z ** 0
        if n == 1: return z - 0.5
        if n == 2: return z*(z-1)
    if abs(z) == ctx.inf:
        return z**n
    if z != z:
        return z
    m = n+1
    if z == 0:
        return -2*(ctx.ldexp(1,m)-1)*ctx.bernoulli(m)/m * z**0
    if z == 1:
        return 2*(ctx.ldexp(1,m)-1)*ctx.bernoulli(m)/m * z**0
    if z == 0.5:
        if n % 2:
            return ctx.zero
        # Use exact code for Euler numbers
        if n < 100 or n*ctx.mag(0.46839865*n) < ctx.prec*0.25:
            return ctx.ldexp(ctx._eulernum(n), -n)
    # http://functions.wolfram.com/Polynomials/EulerE2/06/01/02/01/0002/
    def terms():
        t = ctx.one
        k = 0
        w = ctx.ldexp(1,n+2)
        while 1:
            v = n-k+1
            if not (v > 2 and v & 1):
                yield (2-w)*ctx.bernoulli(v)*t
            k += 1
            if k > n:
                break
            t = t*z*(n-k+2)/k
            w *= 0.5
    return ctx.sum_accurately(terms) / m

@defun
def eulernum(ctx, n, exact=False):
    n = int(n)
    if exact:
        return int(ctx._eulernum(n))
    if n < 100:
        return ctx.mpf(ctx._eulernum(n))
    if n % 2:
        return ctx.zero
    return ctx.ldexp(ctx.eulerpoly(n,0.5), n)

# TODO: this should be implemented low-level
def polylog_series(ctx, s, z):
    tol = +ctx.eps
    l = ctx.zero
    k = 1
    zk = z
    while 1:
        term = zk / k**s
        l += term
        if abs(term) < tol:
            break
        zk *= z
        k += 1
    return l

def polylog_continuation(ctx, n, z):
    if n < 0:
        return z*0
    twopij = 2j * ctx.pi
    a = -twopij**n/ctx.fac(n) * ctx.bernpoly(n, ctx.ln(z)/twopij)
    if ctx._is_real_type(z) and z < 0:
        a = ctx._re(a)
    if ctx._im(z) < 0 or (ctx._im(z) == 0 and ctx._re(z) >= 1):
        a -= twopij*ctx.ln(z)**(n-1)/ctx.fac(n-1)
    return a

def polylog_unitcircle(ctx, n, z):
    tol = +ctx.eps
    if n > 1:
        l = ctx.zero
        logz = ctx.ln(z)
        logmz = ctx.one
        m = 0
        while 1:
            if (n-m) != 1:
                term = ctx.zeta(n-m) * logmz / ctx.fac(m)
                if term and abs(term) < tol:
                    break
                l += term
            logmz *= logz
            m += 1
        l += ctx.ln(z)**(n-1)/ctx.fac(n-1)*(ctx.harmonic(n-1)-ctx.ln(-ctx.ln(z)))
    elif n < 1:  # else
        l = ctx.fac(-n)*(-ctx.ln(z))**(n-1)
        logz = ctx.ln(z)
        logkz = ctx.one
        k = 0
        while 1:
            b = ctx.bernoulli(k-n+1)
            if b:
                term = b*logkz/(ctx.fac(k)*(k-n+1))
                if abs(term) < tol:
                    break
                l -= term
            logkz *= logz
            k += 1
    else:
        raise ValueError
    if ctx._is_real_type(z) and z < 0:
        l = ctx._re(l)
    return l

def polylog_general(ctx, s, z):
    v = ctx.zero
    u = ctx.ln(z)
    if not abs(u) < 5: # theoretically |u| < 2*pi
        raise NotImplementedError("polylog for arbitrary s and z")
    t = 1
    k = 0
    while 1:
        term = ctx.zeta(s-k) * t
        if abs(term) < ctx.eps:
            break
        v += term
        k += 1
        t *= u
        t /= k
    return ctx.gamma(1-s)*(-u)**(s-1) + v

@defun_wrapped
def polylog(ctx, s, z):
    s = ctx.convert(s)
    z = ctx.convert(z)
    if z == 1:
        return ctx.zeta(s)
    if z == -1:
        return -ctx.altzeta(s)
    if s == 0:
        return z/(1-z)
    if s == 1:
        return -ctx.ln(1-z)
    if s == -1:
        return z/(1-z)**2
    if abs(z) <= 0.75 or (not ctx.isint(s) and abs(z) < 0.9):
        return polylog_series(ctx, s, z)
    if abs(z) >= 1.4 and ctx.isint(s):
        return (-1)**(s+1)*polylog_series(ctx, s, 1/z) + polylog_continuation(ctx, s, z)
    if ctx.isint(s):
        return polylog_unitcircle(ctx, int(s), z)
    return polylog_general(ctx, s, z)

    #raise NotImplementedError("polylog for arbitrary s and z")
    # This could perhaps be used in some cases
    #from quadrature import quad
    #return quad(lambda t: t**(s-1)/(exp(t)/z-1),[0,inf])/gamma(s)

@defun_wrapped
def clsin(ctx, s, z, pi=False):
    if ctx.isint(s) and s < 0 and int(s) % 2 == 1:
        return z*0
    if pi:
        a = ctx.expjpi(z)
    else:
        a = ctx.expj(z)
    if ctx._is_real_type(z) and ctx._is_real_type(s):
        return ctx.im(ctx.polylog(s,a))
    b = 1/a
    return (-0.5j)*(ctx.polylog(s,a) - ctx.polylog(s,b))

@defun_wrapped
def clcos(ctx, s, z, pi=False):
    if ctx.isint(s) and s < 0 and int(s) % 2 == 0:
        return z*0
    if pi:
        a = ctx.expjpi(z)
    else:
        a = ctx.expj(z)
    if ctx._is_real_type(z) and ctx._is_real_type(s):
        return ctx.re(ctx.polylog(s,a))
    b = 1/a
    return 0.5*(ctx.polylog(s,a) + ctx.polylog(s,b))

@defun
def altzeta(ctx, s, **kwargs):
    try:
        return ctx._altzeta(s, **kwargs)
    except NotImplementedError:
        return ctx._altzeta_generic(s)

@defun_wrapped
def _altzeta_generic(ctx, s):
    if s == 1:
        return ctx.ln2 + 0*s
    return -ctx.powm1(2, 1-s) * ctx.zeta(s)

@defun
def zeta(ctx, s, a=1, derivative=0, method=None, **kwargs):
    d = int(derivative)
    if a == 1 and not (d or method):
        try:
            return ctx._zeta(s, **kwargs)
        except NotImplementedError:
            pass
    s = ctx.convert(s)
    prec = ctx.prec
    method = kwargs.get('method')
    verbose = kwargs.get('verbose')
    if a == 1 and method != 'euler-maclaurin':
        im = abs(ctx._im(s))
        re = abs(ctx._re(s))
        #if (im < prec or method == 'borwein') and not derivative:
        #    try:
        #        if verbose:
        #            print "zeta: Attempting to use the Borwein algorithm"
        #        return ctx._zeta(s, **kwargs)
        #    except NotImplementedError:
        #        if verbose:
        #            print "zeta: Could not use the Borwein algorithm"
        #        pass
        if abs(im) > 60*prec and 10*re < prec and derivative <= 4 or \
            method == 'riemann-siegel':
            try:   #  py2.4 compatible try block
                try:
                    if verbose:
                        print "zeta: Attempting to use the Riemann-Siegel algorithm"
                    return ctx.rs_zeta(s, derivative, **kwargs)
                except NotImplementedError:
                    if verbose:
                        print "zeta: Could not use the Riemann-Siegel algorithm"
                    pass
            finally:
                ctx.prec = prec
    if s == 1:
        return ctx.inf
    abss = abs(s)
    if abss == ctx.inf:
        if ctx.re(s) == ctx.inf:
            if d == 0:
                return ctx.one
            return ctx.zero
        return s*0
    elif ctx.isnan(abss):
        return 1/s
    if ctx.re(s) > 2*ctx.prec and a == 1 and not derivative:
        return ctx.one + ctx.power(2, -s)
    if verbose:
        print "zeta: Using the Euler-Maclaurin algorithm"
    prec = ctx.prec
    try:
        ctx.prec += 10
        v = ctx._hurwitz(s, a, d)
    finally:
        ctx.prec = prec
    return +v

@defun
def _hurwitz(ctx, s, a=1, d=0):
    # We strongly want to special-case rational a
    a, atype = ctx._convert_param(a)
    prec = ctx.prec
    # TODO: implement reflection for derivatives
    res = ctx.re(s)
    negs = -s
    try:
        if res < 0 and not d:
            # Integer reflection formula
            if ctx.isnpint(s):
                n = int(res)
                if n <= 0:
                    return ctx.bernpoly(1-n, a) / (n-1)
            t = 1-s
            # We now require a to be standardized
            v = 0
            shift = 0
            b = a
            while ctx.re(b) > 1:
                b -= 1
                v -= b**negs
                shift -= 1
            while ctx.re(b) <= 0:
                v += b**negs
                b += 1
                shift += 1
            # Rational reflection formula
            if atype == 'Q' or atype == 'Z':
                try:
                    p, q = a
                except:
                    assert a == int(a)
                    p = int(a)
                    q = 1
                p += shift*q
                assert 1 <= p <= q
                g = ctx.fsum(ctx.cospi(t/2-2*k*b)*ctx._hurwitz(t,(k,q)) \
                    for k in range(1,q+1))
                g *= 2*ctx.gamma(t)/(2*ctx.pi*q)**t
                v += g
                return v
            # General reflection formula
            else:
                C1 = ctx.cospi(t/2)
                C2 = ctx.sinpi(t/2)
                # Clausen functions; could maybe use polylog directly
                if C1: C1 *= ctx.clcos(t, 2*a, pi=True)
                if C2: C2 *= ctx.clsin(t, 2*a, pi=True)
                v += 2*ctx.gamma(t)/(2*ctx.pi)**t*(C1+C2)
                return v
    except NotImplementedError:
        pass
    a = ctx.convert(a)
    tol = -prec
    # Estimate number of terms for Euler-Maclaurin summation; could be improved
    M1 = 0
    M2 = prec // 3
    N = M2
    lsum = 0
    # This speeds up the recurrence for derivatives
    if ctx.isint(s):
        s = int(ctx._re(s))
    s1 = s-1
    while 1:
        # Truncated L-series
        l = ctx._zetasum(s, M1+a, M2-M1-1, [d])[0][0]
        #if d:
        #    l = ctx.fsum((-ctx.ln(n+a))**d * (n+a)**negs for n in range(M1,M2))
        #else:
        #    l = ctx.fsum((n+a)**negs for n in range(M1,M2))
        lsum += l
        M2a = M2+a
        logM2a = ctx.ln(M2a)
        logM2ad = logM2a**d
        logs = [logM2ad]
        logr = 1/logM2a
        rM2a = 1/M2a
        M2as = rM2a**s
        if d:
            tailsum = ctx.gammainc(d+1, s1*logM2a) / s1**(d+1)
        else:
            tailsum = 1/((s1)*(M2a)**s1)
        tailsum += 0.5 * logM2ad * M2as
        U = [1]
        r = M2as
        fact = 2
        for j in range(1, N+1):
            # TODO: the following could perhaps be tidied a bit
            j2 = 2*j
            if j == 1:
                upds = [1]
            else:
                upds = [j2-2, j2-1]
            for m in upds:
                D = min(m,d+1)
                if m <= d:
                    logs.append(logs[-1] * logr)
                Un = [0]*(D+1)
                for i in xrange(D): Un[i] = (1-m-s)*U[i]
                for i in xrange(1,D+1): Un[i] += (d-(i-1))*U[i-1]
                U = Un
                r *= rM2a
            t = ctx.fdot(U, logs) * r * ctx.bernoulli(j2)/(-fact)
            tailsum += t
            if ctx.mag(t) < tol:
                return lsum + (-1)**d * tailsum
            fact *= (j2+1)*(j2+2)
        M1, M2 = M2, M2*2

@defun
def _zetasum(ctx, s, a, n, derivatives=[0], reflect=False):
    """
    Returns [xd0,xd1,...,xdr], [yd0,yd1,...ydr] where

    xdk = D^k     ( 1/a^s     +  1/(a+1)^s      +  ...  +  1/(a+n)^s     )
    ydk = D^k conj( 1/a^(1-s) +  1/(a+1)^(1-s)  +  ...  +  1/(a+n)^(1-s) )

    D^k = kth derivative with respect to s, k ranges over the given list of
    derivatives (which should consist of either a single element
    or a range 0,1,...r). If reflect=False, the ydks are not computed.
    """
    try:
        return ctx._zetasum_fast(s, a, n, derivatives, reflect)
    except NotImplementedError:
        pass
    negs = ctx.fneg(s, exact=True)
    have_derivatives = derivatives != [0]
    have_one_derivative = len(derivatives) == 1
    if not reflect:
        if not have_derivatives:
            return [ctx.fsum((a+k)**negs for k in xrange(n+1))], []
        if have_one_derivative:
            d = derivatives[0]
            x = ctx.fsum(ctx.ln(a+k)**d * (a+k)**negs for k in xrange(n+1))
            return [(-1)**d * x], []
    maxd = max(derivatives)
    if not have_one_derivative:
        derivatives = range(maxd+1)
    xs = [ctx.zero for d in derivatives]
    if reflect:
        ys = [ctx.zero for d in derivatives]
    else:
        ys = []
    for k in xrange(n+1):
        w = a + k
        xterm = w ** negs
        if reflect:
            yterm = ctx.conj(ctx.one / (w * xterm))
        if have_derivatives:
            logw = -ctx.ln(w)
            if have_one_derivative:
                logw = logw ** maxd
                xs[0] += xterm * logw
                if reflect:
                    ys[0] += yterm * logw
            else:
                t = ctx.one
                for d in derivatives:
                    xs[d] += xterm * t
                    if reflect:
                        ys[d] += yterm * t
                    t *= logw
        else:
            xs[0] += xterm
            if reflect:
                ys[0] += yterm
    return xs, ys

@defun
def dirichlet(ctx, s, chi=[1], derivative=0):
    s = ctx.convert(s)
    q = len(chi)
    d = int(derivative)
    if d > 2:
        raise NotImplementedError("arbitrary order derivatives")
    prec = ctx.prec
    try:
        ctx.prec += 10
        if s == 1:
            have_pole = True
            for x in chi:
                if x and x != 1:
                    have_pole = False
                    h = +ctx.eps
                    ctx.prec *= 2*(d+1)
                    s += h
            if have_pole:
                return +ctx.inf
        z = ctx.zero
        for p in range(1,q+1):
            if chi[p%q]:
                if d == 1:
                    z += chi[p%q] * (ctx.zeta(s, (p,q), 1) - \
                        ctx.zeta(s, (p,q))*ctx.log(q))
                else:
                    z += chi[p%q] * ctx.zeta(s, (p,q))
        z /= q**s
    finally:
        ctx.prec = prec
    return +z
