from functions import defun, defun_wrapped

def _check_need_perturb(ctx, terms, prec, discard_known_zeros):
    perturb = recompute = False
    extraprec = 0
    discard = []
    for term_index, term in enumerate(terms):
        w_s, c_s, alpha_s, beta_s, a_s, b_s, z = term
        have_singular_nongamma_weight = False
        # Avoid division by zero in leading factors (TODO:
        # also check for near division by zero?)
        for k, w in enumerate(w_s):
            if not w:
                if ctx.re(c_s[k]) <= 0 and c_s[k]:
                    perturb = recompute = True
                    have_singular_nongamma_weight = True
        pole_count = [0, 0, 0]
        # Check for gamma and series poles and near-poles
        for data_index, data in enumerate([alpha_s, beta_s, b_s]):
            for i, x in enumerate(data):
                n, d = ctx.nint_distance(x)
                # Poles
                if n > 0:
                    continue
                if d == ctx.ninf:
                    # OK if we have a polynomial
                    # ------------------------------
                    ok = False
                    if data_index == 2:
                        for u in a_s:
                            if ctx.isnpint(u) and u >= int(n):
                                ok = True
                                break
                    if ok:
                        continue
                    pole_count[data_index] += 1
                    # ------------------------------
                    #perturb = recompute = True
                    #return perturb, recompute, extraprec
                elif d < -4:
                    extraprec += -d
                    recompute = True
        if discard_known_zeros and pole_count[1] > pole_count[0] + pole_count[2] \
            and not have_singular_nongamma_weight:
            discard.append(term_index)
        elif sum(pole_count):
            perturb = recompute = True
    return perturb, recompute, extraprec, discard

_hypercomb_msg = """
hypercomb() failed to converge to the requested %i bits of accuracy
using a working precision of %i bits. The function value may be zero or
infinite; try passing zeroprec=N or infprec=M to bound finite values between
2^(-N) and 2^M. Otherwise try a higher maxprec or maxterms.
"""

@defun
def hypercomb(ctx, function, params=[], discard_known_zeros=True, **kwargs):
    orig = ctx.prec
    sumvalue = ctx.zero
    dist = ctx.nint_distance
    ninf = ctx.ninf
    orig_params = params[:]
    verbose = kwargs.get('verbose', False)
    maxprec = kwargs.get('maxprec', ctx._default_hyper_maxprec(orig))
    kwargs['maxprec'] = maxprec   # For calls to hypsum
    zeroprec = kwargs.get('zeroprec')
    infprec = kwargs.get('infprec')
    perturbed_reference_value = None
    hextra = 0
    try:
        while 1:
            ctx.prec += 10
            if ctx.prec > maxprec:
                raise ValueError(_hypercomb_msg % (orig, ctx.prec))
            orig2 = ctx.prec
            params = orig_params[:]
            terms = function(*params)
            if verbose:
                print
                print "ENTERING hypercomb main loop"
                print "prec =", ctx.prec
                print "hextra", hextra
            perturb, recompute, extraprec, discard = \
                _check_need_perturb(ctx, terms, orig, discard_known_zeros)
            ctx.prec += extraprec
            if perturb:
                if "hmag" in kwargs:
                    hmag = kwargs["hmag"]
                elif ctx._fixed_precision:
                    hmag = int(ctx.prec*0.3)
                else:
                    hmag = orig + 10 + hextra
                h = ctx.ldexp(ctx.one, -hmag)
                ctx.prec = orig2 + 10 + hmag + 10
                for k in range(len(params)):
                    params[k] += h
                    # Heuristically ensure that the perturbations
                    # are "independent" so that two perturbations
                    # don't accidentally cancel each other out
                    # in a subtraction.
                    h += h/(k+1)
            if recompute:
                terms = function(*params)
            if discard_known_zeros:
                terms = [term for (i, term) in enumerate(terms) if i not in discard]
            if not terms:
                return ctx.zero
            evaluated_terms = []
            for term_index, term_data in enumerate(terms):
                w_s, c_s, alpha_s, beta_s, a_s, b_s, z = term_data
                if verbose:
                    print
                    print "  Evaluating term %i/%i : %iF%i" % \
                        (term_index+1, len(terms), len(a_s), len(b_s))
                    print "    powers", ctx.nstr(w_s), ctx.nstr(c_s)
                    print "    gamma", ctx.nstr(alpha_s), ctx.nstr(beta_s)
                    print "    hyper", ctx.nstr(a_s), ctx.nstr(b_s)
                    print "    z", ctx.nstr(z)
                v = ctx.hyper(a_s, b_s, z, **kwargs)
                for a in alpha_s: v *= ctx.gamma(a)
                for b in beta_s: v /= ctx.gamma(b)
                for w, c in zip(w_s, c_s): v *= ctx.power(w, c)
                if verbose:
                    print "    Value:", v
                evaluated_terms.append(v)

            if len(terms) == 1 and (not perturb):
                sumvalue = evaluated_terms[0]
                break

            if ctx._fixed_precision:
                sumvalue = ctx.fsum(evaluated_terms)
                break

            sumvalue = ctx.fsum(evaluated_terms)
            term_magnitudes = [ctx.mag(x) for x in evaluated_terms]
            max_magnitude = max(term_magnitudes)
            sum_magnitude = ctx.mag(sumvalue)
            cancellation = max_magnitude - sum_magnitude
            if verbose:
                print
                print "  Cancellation:", cancellation, "bits"
                print "  Increased precision:", ctx.prec - orig, "bits"

            precision_ok = cancellation < ctx.prec - orig

            if zeroprec is None:
                zero_ok = False
            else:
                zero_ok = max_magnitude - ctx.prec < -zeroprec
            if infprec is None:
                inf_ok = False
            else:
                inf_ok = max_magnitude > infprec

            if precision_ok and (not perturb) or ctx.isnan(cancellation):
                break
            elif precision_ok:
                if perturbed_reference_value is None:
                    hextra += 20
                    perturbed_reference_value = sumvalue
                    continue
                elif ctx.mag(sumvalue - perturbed_reference_value) <= \
                        ctx.mag(sumvalue) - orig:
                    break
                elif zero_ok:
                    sumvalue = ctx.zero
                    break
                elif inf_ok:
                    sumvalue = ctx.inf
                    break
                elif 'hmag' in kwargs:
                    break
                else:
                    hextra *= 2
                    perturbed_reference_value = sumvalue
            # Increase precision
            else:
                increment = min(max(cancellation, orig//2), max(extraprec,orig))
                ctx.prec += increment
                if verbose:
                    print "  Must start over with increased precision"
                continue
    finally:
        ctx.prec = orig
    return +sumvalue

@defun
def hyper(ctx, a_s, b_s, z, **kwargs):
    """
    Hypergeometric function, general case.
    """
    z = ctx.convert(z)
    p = len(a_s)
    q = len(b_s)
    a_s = map(ctx._convert_param, a_s)
    b_s = map(ctx._convert_param, b_s)
    # Reduce degree by eliminating common parameters
    if kwargs.get('eliminate', True):
        i = 0
        while i < q and a_s:
            b = b_s[i]
            if b in a_s:
                a_s.remove(b)
                b_s.remove(b)
                p -= 1
                q -= 1
            else:
                i += 1
    # Handle special cases
    if p == 0:
        if   q == 1: return ctx._hyp0f1(b_s, z, **kwargs)
        elif q == 0: return ctx.exp(z)
    elif p == 1:
        if   q == 1: return ctx._hyp1f1(a_s, b_s, z, **kwargs)
        elif q == 2: return ctx._hyp1f2(a_s, b_s, z, **kwargs)
        elif q == 0: return ctx._hyp1f0(a_s[0][0], z)
    elif p == 2:
        if   q == 1: return ctx._hyp2f1(a_s, b_s, z, **kwargs)
        elif q == 2: return ctx._hyp2f2(a_s, b_s, z, **kwargs)
        elif q == 3: return ctx._hyp2f3(a_s, b_s, z, **kwargs)
        elif q == 0: return ctx._hyp2f0(a_s, b_s, z, **kwargs)
    elif p == q+1:
        return ctx._hypq1fq(p, q, a_s, b_s, z, **kwargs)
    elif p > q+1 and not kwargs.get('force_series'):
        return ctx._hyp_borel(p, q, a_s, b_s, z, **kwargs)
    coeffs, types = zip(*(a_s+b_s))
    return ctx.hypsum(p, q, types, coeffs, z, **kwargs)

@defun
def hyp0f1(ctx,b,z,**kwargs):
    return ctx.hyper([],[b],z,**kwargs)

@defun
def hyp1f1(ctx,a,b,z,**kwargs):
    return ctx.hyper([a],[b],z,**kwargs)

@defun
def hyp1f2(ctx,a1,b1,b2,z,**kwargs):
    return ctx.hyper([a1],[b1,b2],z,**kwargs)

@defun
def hyp2f1(ctx,a,b,c,z,**kwargs):
    return ctx.hyper([a,b],[c],z,**kwargs)

@defun
def hyp2f2(ctx,a1,a2,b1,b2,z,**kwargs):
    return ctx.hyper([a1,a2],[b1,b2],z,**kwargs)

@defun
def hyp2f3(ctx,a1,a2,b1,b2,b3,z,**kwargs):
    return ctx.hyper([a1,a2],[b1,b2,b3],z,**kwargs)

@defun
def hyp2f0(ctx,a,b,z,**kwargs):
    return ctx.hyper([a,b],[],z,**kwargs)

@defun
def hyp3f2(ctx,a1,a2,a3,b1,b2,z,**kwargs):
    return ctx.hyper([a1,a2,a3],[b1,b2],z,**kwargs)

@defun_wrapped
def _hyp1f0(ctx, a, z):
    return (1-z) ** (-a)

@defun
def _hyp0f1(ctx, b_s, z, **kwargs):
    (b, btype), = b_s
    if z:
        magz = ctx.mag(z)
    else:
        magz = 0
    if magz >= 8 and not kwargs.get('force_series'):
        try:
            # http://functions.wolfram.com/HypergeometricFunctions/
            # Hypergeometric0F1/06/02/03/0004/
            # We don't need hypercomb because the only possible singularity
            # occurs when the value is undefined. However, we should perhaps
            # still check for cancellation...
            # TODO: handle the all-real case more efficiently!
            # TODO: figure out how much precision is needed (exponential growth)
            orig = ctx.prec
            try:
                ctx.prec += 12 + magz//2
                w = ctx.sqrt(-z)
                jw = ctx.j*w
                u = 1/(4*jw)
                c = ctx.mpq_1_2 - b
                E = ctx.exp(2*jw)
                H1 = (-jw)**c/E*ctx.hyp2f0(b-ctx.mpq_1_2, ctx.mpq_3_2-b, -u,
                    force_series=True)
                H2 = (jw)**c*E*ctx.hyp2f0(b-ctx.mpq_1_2, ctx.mpq_3_2-b, u,
                    force_series=True)
                v = ctx.gamma(b)/(2*ctx.sqrt(ctx.pi))*(H1 + H2)
            finally:
                ctx.prec = orig
            if ctx._is_real_type(b) and ctx._is_real_type(z):
                v = ctx._re(v)
            return +v
        except ctx.NoConvergence:
            pass
    return ctx.hypsum(0, 1, (btype,), [b], z, **kwargs)

@defun
def _hyp1f1(ctx, a_s, b_s, z, **kwargs):
    (a, atype), = a_s
    (b, btype), = b_s
    if not z:
        return ctx.one+z
    magz = ctx.mag(z)
    if magz >= 7 and not (ctx.isint(a) and ctx.re(a) <= 0):
        if ctx.isinf(z):
            if ctx.sign(a) == ctx.sign(b) == ctx.sign(z) == 1:
                return ctx.inf
            return ctx.nan * z
        try:
            try:
                ctx.prec += magz
                sector = ctx._im(z) < 0 and ctx._re(z) <= 0
                def h(a,b):
                    if sector:
                        E = ctx.expjpi(ctx.fneg(a, exact=True))
                    else:
                        E = ctx.expjpi(a)
                    rz = 1/z
                    T1 = ([E,z], [1,-a], [b], [b-a], [a, 1+a-b], [], -rz)
                    T2 = ([ctx.exp(z),z], [1,a-b], [b], [a], [b-a, 1-a], [], rz)
                    return T1, T2
                v = ctx.hypercomb(h, [a,b], force_series=True)
                if ctx._is_real_type(a) and ctx._is_real_type(b) and ctx._is_real_type(z):
                    v = ctx._re(v)
                return +v
            except ctx.NoConvergence:
                pass
        finally:
            ctx.prec -= magz
    v = ctx.hypsum(1, 1, (atype, btype), [a, b], z, **kwargs)
    return v

def _hyp2f1_gosper(ctx,a,b,c,z,**kwargs):
    # Use Gosper's recurrence
    # See http://www.math.utexas.edu/pipermail/maxima/2006/000126.html
    _a,_b,_c,_z = a, b, c, z
    orig = ctx.prec
    maxprec = kwargs.get('maxprec', 100*orig)
    extra = 10
    while 1:
        ctx.prec = orig + extra
        #a = ctx.convert(_a)
        #b = ctx.convert(_b)
        #c = ctx.convert(_c)
        z = ctx.convert(_z)
        d = ctx.mpf(0)
        e = ctx.mpf(1)
        f = ctx.mpf(0)
        k = 0
        # Common subexpression elimination, unfortunately making
        # things a bit unreadable. The formula is quite messy to begin
        # with, though...
        abz = a*b*z
        ch = c * ctx.mpq_1_2
        c1h = (c+1) * ctx.mpq_1_2
        nz = 1-z
        g = z/nz
        abg = a*b*g
        cba = c-b-a
        z2 = z-2
        tol = -ctx.prec - 10
        nstr = ctx.nstr
        nprint = ctx.nprint
        mag = ctx.mag
        maxmag = ctx.ninf
        while 1:
            kch = k+ch
            kakbz = (k+a)*(k+b)*z / (4*(k+1)*kch*(k+c1h))
            d1 = kakbz*(e-(k+cba)*d*g)
            e1 = kakbz*(d*abg+(k+c)*e)
            ft = d*(k*(cba*z+k*z2-c)-abz)/(2*kch*nz)
            f1 = f + e - ft
            maxmag = max(maxmag, mag(f1))
            if mag(f1-f) < tol:
                break
            d, e, f = d1, e1, f1
            k += 1
        cancellation = maxmag - mag(f1)
        if cancellation < extra:
            break
        else:
            extra += cancellation
            if extra > maxprec:
                raise ctx.NoConvergence
    return f1

@defun
def _hyp2f1(ctx, a_s, b_s, z, **kwargs):
    (a, atype), (b, btype) = a_s
    (c, ctype), = b_s
    if z == 1:
        # TODO: the following logic can be simplified
        convergent = ctx.re(c-a-b) > 0
        finite = (ctx.isint(a) and a <= 0) or (ctx.isint(b) and b <= 0)
        zerodiv = ctx.isint(c) and c <= 0 and not \
            ((ctx.isint(a) and c <= a <= 0) or (ctx.isint(b) and c <= b <= 0))
        #print "bz", a, b, c, z, convergent, finite, zerodiv
        # Gauss's theorem gives the value if convergent
        if (convergent or finite) and not zerodiv:
            return ctx.gammaprod([c, c-a-b], [c-a, c-b], _infsign=True)
        # Otherwise, there is a pole and we take the
        # sign to be that when approaching from below
        # XXX: this evaluation is not necessarily correct in all cases
        return ctx.hyp2f1(a,b,c,1-ctx.eps*2) * ctx.inf

    # Equal to 1 (first term), unless there is a subsequent
    # division by zero
    if not z:
        # Division by zero but power of z is higher than
        # first order so cancels
        if c or a == 0 or b == 0:
            return 1+z
        # Indeterminate
        return ctx.nan

    # Hit zero denominator unless numerator goes to 0 first
    if ctx.isint(c) and c <= 0:
        if (ctx.isint(a) and c <= a <= 0) or \
           (ctx.isint(b) and c <= b <= 0):
            pass
        else:
            # Pole in series
            return ctx.inf

    absz = abs(z)

    # Fast case: standard series converges rapidly,
    # possibly in finitely many terms
    if absz <= 0.8 or (ctx.isint(a) and a <= 0 and a >= -1000) or \
                      (ctx.isint(b) and b <= 0 and b >= -1000):
        return ctx.hypsum(2, 1, (atype, btype, ctype), [a, b, c], z, **kwargs)

    orig = ctx.prec
    try:
        ctx.prec += 10

        # Use 1/z transformation
        if absz >= 1.3:
            def h(a,b):
                t = ctx.mpq_1-c; ab = a-b; rz = 1/z
                T1 = ([-z],[-a], [c,-ab],[b,c-a], [a,t+a],[ctx.mpq_1+ab],  rz)
                T2 = ([-z],[-b], [c,ab],[a,c-b], [b,t+b],[ctx.mpq_1-ab],  rz)
                return T1, T2
            v = ctx.hypercomb(h, [a,b], **kwargs)

        # Use 1-z transformation
        elif abs(1-z) <= 0.75:
            def h(a,b):
                t = c-a-b; ca = c-a; cb = c-b; rz = 1-z
                T1 = [], [], [c,t], [ca,cb], [a,b], [1-t], rz
                T2 = [rz], [t], [c,a+b-c], [a,b], [ca,cb], [1+t], rz
                return T1, T2
            v = ctx.hypercomb(h, [a,b], **kwargs)

        # Use z/(z-1) transformation
        elif abs(z/(z-1)) <= 0.75:
            v = ctx.hyp2f1(a, c-b, c, z/(z-1)) / (1-z)**a

        # Remaining part of unit circle
        else:
            v = _hyp2f1_gosper(ctx,a,b,c,z,**kwargs)

    finally:
        ctx.prec = orig
    return +v

@defun
def _hypq1fq(ctx, p, q, a_s, b_s, z, **kwargs):
    r"""
    Evaluates 3F2, 4F3, 5F4, ...
    """
    a_s, a_types = zip(*a_s)
    b_s, b_types = zip(*b_s)
    a_s = list(a_s)
    b_s = list(b_s)
    absz = abs(z)
    ispoly = False
    for a in a_s:
        if ctx.isint(a) and a <= 0:
            ispoly = True
            break
    # Direct summation
    if absz < 1 or ispoly:
        try:
            return ctx.hypsum(p, q, a_types+b_types, a_s+b_s, z, **kwargs)
        except ctx.NoConvergence:
            if absz > 1.1 or ispoly:
                raise
    # Use expansion at |z-1| -> 0.
    # Reference: Wolfgang Buhring, "Generalized Hypergeometric Functions at
    #   Unit Argument", Proc. Amer. Math. Soc., Vol. 114, No. 1 (Jan. 1992),
    #   pp.145-153
    # The current implementation has several problems:
    # 1. We only implement it for 3F2. The expansion coefficients are
    #    given by extremely messy nested sums in the higher degree cases
    #    (see reference). Is efficient sequential generation of the coefficients
    #    possible in the > 3F2 case?
    # 2. Although the series converges, it may do so slowly, so we need
    #    convergence acceleration. The acceleration implemented by
    #    nsum does not always help, so results returned are sometimes
    #    inaccurate! Can we do better?
    # 3. We should check conditions for convergence, and possibly
    #    do a better job of cancelling out gamma poles if possible.
    if z == 1:
        # XXX: should also check for division by zero in the
        # denominator of the series (cf. hyp2f1)
        S = ctx.re(sum(b_s)-sum(a_s))
        if S <= 0:
            #return ctx.hyper(a_s, b_s, 1-ctx.eps*2, **kwargs) * ctx.inf
            return ctx.hyper(a_s, b_s, 0.9, **kwargs) * ctx.inf
    if (p,q) == (3,2) and abs(z-1) < 0.05:   # and kwargs.get('sum1')
        #print "Using alternate summation (experimental)"
        a1,a2,a3 = a_s
        b1,b2 = b_s
        u = b1+b2-a3
        initial = ctx.gammaprod([b2-a3,b1-a3,a1,a2],[b2-a3,b1-a3,1,u])
        def term(k, _cache={0:initial}):
            u = b1+b2-a3+k
            if k in _cache:
                t = _cache[k]
            else:
                t = _cache[k-1]
                t *= (b1+k-a3-1)*(b2+k-a3-1)
                t /= k*(u-1)
                _cache[k] = t
            return t * ctx.hyp2f1(a1,a2,u,z)
        try:
            S = ctx.nsum(term, [0,ctx.inf], verbose=kwargs.get('verbose'),
                strict=kwargs.get('strict', True))
            return S * ctx.gammaprod([b1,b2],[a1,a2,a3])
        except ctx.NoConvergence:
            pass
    # Try to use convergence acceleration on and close to the unit circle.
    # Problem: the convergence acceleration degenerates as |z-1| -> 0,
    # except for special cases. Everywhere else, the Shanks transformation
    # is very efficient.
    if absz < 1.1 and ctx._re(z) <= 1:
        def term(k, _cache={0:ctx.one}):
            k = int(k)
            if k in _cache:
                return _cache[k]
            t = _cache[k-1]
            m = k-1
            for j in xrange(p): t *= (a_s[j]+m)
            for j in xrange(q): t /= (b_s[j]+m)
            t *= z
            t /= k
            _cache[k] = t
            return t
        return ctx.nsum(term, [0,ctx.inf], verbose=kwargs.get('verbose'),
            strict=kwargs.get('strict', True))
    # Use 1/z transformation
    # http://functions.wolfram.com/HypergeometricFunctions/
    #   HypergeometricPFQ/06/01/05/02/0004/
    def h(*args):
        a_s = list(args[:p])
        b_s = list(args[p:])
        Ts = []
        recz = ctx.one/z
        negz = ctx.fneg(z, exact=True)
        for k in range(q+1):
            ak = a_s[k]
            C = [negz]
            Cp = [-ak]
            Gn = b_s + [ak] + [a_s[j]-ak for j in range(q+1) if j != k]
            Gd = a_s + [b_s[j]-ak for j in range(q)]
            Fn = [ak] + [ak-b_s[j]+1 for j in range(q)]
            Fd = [1-a_s[j]+ak for j in range(q+1) if j != k]
            Ts.append((C, Cp, Gn, Gd, Fn, Fd, recz))
        return Ts
    return ctx.hypercomb(h, a_s+b_s, **kwargs)

@defun
def _hyp_borel(ctx, p, q, a_s, b_s, z, **kwargs):
    if a_s:
        a_s, a_types = zip(*a_s)
        a_s = list(a_s)
    else:
        a_s, a_types = [], ()
    if b_s:
        b_s, b_types = zip(*b_s)
        b_s = list(b_s)
    else:
        b_s, b_types = [], ()
    kwargs['maxterms'] = kwargs.get('maxterms', ctx.prec)
    try:
        return ctx.hypsum(p, q, a_types+b_types, a_s+b_s, z, **kwargs)
    except ctx.NoConvergence:
        pass
    prec = ctx.prec
    try:
        tol = kwargs.get('asymp_tol', ctx.eps/4)
        ctx.prec += 10
        # hypsum is has a conservative tolerance. So we try again:
        def term(k, cache={0:ctx.one}):
            if k in cache:
                return cache[k]
            t = term(k-1)
            for a in a_s: t *= (a+(k-1))
            for b in b_s: t /= (b+(k-1))
            t *= z
            t /= k
            cache[k] = t
            return t
        s = ctx.one
        for k in xrange(1, ctx.prec):
            t = term(k)
            s += t
            if abs(t) <= tol:
                return s
    finally:
        ctx.prec = prec
    if p <= q+3:
        contour = kwargs.get('contour')
        if not contour:
            if ctx.arg(z) < 0.25:
                u = z / max(1, abs(z))
                if ctx.arg(z) >= 0:
                    contour = [0, 2j, (2j+2)/u, 2/u, ctx.inf]
                else:
                    contour = [0, -2j, (-2j+2)/u, 2/u, ctx.inf]
                #contour = [0, 2j/z, 2/z, ctx.inf]
                #contour = [0, 2j, 2/z, ctx.inf]
                #contour = [0, 2j, ctx.inf]
            else:
                contour = [0, ctx.inf]
        quad_kwargs = kwargs.get('quad_kwargs', {})
        def g(t):
            return ctx.exp(-t)*ctx.hyper(a_s, b_s+[1], t*z)
        I, err = ctx.quad(g, contour, error=True, **quad_kwargs)
        if err <= abs(I)*ctx.eps*8:
            return I
    raise ctx.NoConvergence


@defun
def _hyp2f2(ctx, a_s, b_s, z, **kwargs):
    (a1, a1type), (a2, a2type) = a_s
    (b1, b1type), (b2, b2type) = b_s

    absz = abs(z)
    magz = ctx.mag(z)
    orig = ctx.prec

    # Asymptotic expansion is ~ exp(z)
    asymp_extraprec = magz

    # Asymptotic series is in terms of 3F1
    can_use_asymptotic = (not kwargs.get('force_series')) and \
        (ctx.mag(absz) > 3)

    # TODO: much of the following could be shared with 2F3 instead of
    # copypasted
    if can_use_asymptotic:
        #print "using asymp"
        try:
            try:
                ctx.prec += asymp_extraprec
                # http://functions.wolfram.com/HypergeometricFunctions/
                # Hypergeometric2F2/06/02/02/0002/
                def h(a1,a2,b1,b2):
                    X = a1+a2-b1-b2
                    A2 = a1+a2
                    B2 = b1+b2
                    c = {}
                    c[0] = ctx.one
                    c[1] = (A2-1)*X+b1*b2-a1*a2
                    s1 = 0
                    k = 0
                    tprev = 0
                    while 1:
                        if k not in c:
                            uu1 = 1-B2+2*a1+a1**2+2*a2+a2**2-A2*B2+a1*a2+b1*b2+(2*B2-3*(A2+1))*k+2*k**2
                            uu2 = (k-A2+b1-1)*(k-A2+b2-1)*(k-X-2)
                            c[k] = ctx.one/k * (uu1*c[k-1]-uu2*c[k-2])
                        t1 = c[k] * z**(-k)
                        if abs(t1) < 0.1*ctx.eps:
                            #print "Convergence :)"
                            break
                        # Quit if the series doesn't converge quickly enough
                        if k > 5 and abs(tprev) / abs(t1) < 1.5:
                            #print "No convergence :("
                            raise ctx.NoConvergence
                        s1 += t1
                        tprev = t1
                        k += 1
                    S = ctx.exp(z)*s1
                    T1 = [z,S], [X,1], [b1,b2],[a1,a2],[],[],0
                    T2 = [-z],[-a1],[b1,b2,a2-a1],[a2,b1-a1,b2-a1],[a1,a1-b1+1,a1-b2+1],[a1-a2+1],-1/z
                    T3 = [-z],[-a2],[b1,b2,a1-a2],[a1,b1-a2,b2-a2],[a2,a2-b1+1,a2-b2+1],[-a1+a2+1],-1/z
                    return T1, T2, T3
                v = ctx.hypercomb(h, [a1,a2,b1,b2], force_series=True, maxterms=4*ctx.prec)
                if sum(ctx._is_real_type(u) for u in [a1,a2,b1,b2,z]) == 5:
                    v = ctx.re(v)
                return v
            except ctx.NoConvergence:
                pass
        finally:
            ctx.prec = orig

    return ctx.hypsum(2, 2, (a1type, a2type, b1type, b2type), [a1, a2, b1, b2], z, **kwargs)



@defun
def _hyp1f2(ctx, a_s, b_s, z, **kwargs):
    (a1, a1type), = a_s
    (b1, b1type), (b2, b2type) = b_s

    absz = abs(z)
    magz = ctx.mag(z)
    orig = ctx.prec

    # Asymptotic expansion is ~ exp(sqrt(z))
    asymp_extraprec = z and magz//2

    # Asymptotic series is in terms of 3F0
    can_use_asymptotic = (not kwargs.get('force_series')) and \
        (ctx.mag(absz) > 19) and \
        (ctx.sqrt(absz) > 1.5*orig) #and \
        #ctx._hyp_check_convergence([a1, a1-b1+1, a1-b2+1], [],
        #    1/absz, orig+40+asymp_extraprec)

    # TODO: much of the following could be shared with 2F3 instead of
    # copypasted
    if can_use_asymptotic:
        #print "using asymp"
        try:
            try:
                ctx.prec += asymp_extraprec
                # http://functions.wolfram.com/HypergeometricFunctions/
                # Hypergeometric1F2/06/02/03/
                def h(a1,b1,b2):
                    X = ctx.mpq_1_2*(a1-b1-b2+ctx.mpq_1_2)
                    c = {}
                    c[0] = ctx.one
                    c[1] = 2*(ctx.mpq_1_4*(3*a1+b1+b2-2)*(a1-b1-b2)+b1*b2-ctx.mpq_3_16)
                    c[2] = 2*(b1*b2+ctx.mpq_1_4*(a1-b1-b2)*(3*a1+b1+b2-2)-ctx.mpq_3_16)**2+\
                        ctx.mpq_1_16*(-16*(2*a1-3)*b1*b2 + \
                        4*(a1-b1-b2)*(-8*a1**2+11*a1+b1+b2-2)-3)
                    s1 = 0
                    s2 = 0
                    k = 0
                    tprev = 0
                    while 1:
                        if k not in c:
                            uu1 = (3*k**2+(-6*a1+2*b1+2*b2-4)*k + 3*a1**2 - \
                                (b1-b2)**2 - 2*a1*(b1+b2-2) + ctx.mpq_1_4)
                            uu2 = (k-a1+b1-b2-ctx.mpq_1_2)*(k-a1-b1+b2-ctx.mpq_1_2)*\
                                (k-a1+b1+b2-ctx.mpq_5_2)
                            c[k] = ctx.one/(2*k)*(uu1*c[k-1]-uu2*c[k-2])
                        w = c[k] * (-z)**(-0.5*k)
                        t1 = (-ctx.j)**k * ctx.mpf(2)**(-k) * w
                        t2 = ctx.j**k * ctx.mpf(2)**(-k) * w
                        if abs(t1) < 0.1*ctx.eps:
                            #print "Convergence :)"
                            break
                        # Quit if the series doesn't converge quickly enough
                        if k > 5 and abs(tprev) / abs(t1) < 1.5:
                            #print "No convergence :("
                            raise ctx.NoConvergence
                        s1 += t1
                        s2 += t2
                        tprev = t1
                        k += 1
                    S = ctx.expj(ctx.pi*X+2*ctx.sqrt(-z))*s1 + \
                        ctx.expj(-(ctx.pi*X+2*ctx.sqrt(-z)))*s2
                    T1 = [0.5*S, ctx.pi, -z], [1, -0.5, X], [b1, b2], [a1],\
                        [], [], 0
                    T2 = [-z], [-a1], [b1,b2],[b1-a1,b2-a1], \
                        [a1,a1-b1+1,a1-b2+1], [], 1/z
                    return T1, T2
                v = ctx.hypercomb(h, [a1,b1,b2], force_series=True, maxterms=4*ctx.prec)
                if sum(ctx._is_real_type(u) for u in [a1,b1,b2,z]) == 4:
                    v = ctx.re(v)
                return v
            except ctx.NoConvergence:
                pass
        finally:
            ctx.prec = orig

    #print "not using asymp"
    return ctx.hypsum(1, 2, (a1type, b1type, b2type), [a1, b1, b2], z, **kwargs)



@defun
def _hyp2f3(ctx, a_s, b_s, z, **kwargs):
    (a1, a1type), (a2, a2type) = a_s
    (b1, b1type), (b2, b2type), (b3, b3type) = b_s

    absz = abs(z)
    magz = ctx.mag(z)

    # Asymptotic expansion is ~ exp(sqrt(z))
    asymp_extraprec = z and magz//2
    orig = ctx.prec

    # Asymptotic series is in terms of 4F1
    # The square root below empirically provides a plausible criterion
    # for the leading series to converge
    can_use_asymptotic = (not kwargs.get('force_series')) and \
        (ctx.mag(absz) > 19) and (ctx.sqrt(absz) > 1.5*orig)

    if can_use_asymptotic:
        #print "using asymp"
        try:
            try:
                ctx.prec += asymp_extraprec
                # http://functions.wolfram.com/HypergeometricFunctions/
                # Hypergeometric2F3/06/02/03/01/0002/
                def h(a1,a2,b1,b2,b3):
                    X = ctx.mpq_1_2*(a1+a2-b1-b2-b3+ctx.mpq_1_2)
                    A2 = a1+a2
                    B3 = b1+b2+b3
                    A = a1*a2
                    B = b1*b2+b3*b2+b1*b3
                    R = b1*b2*b3
                    c = {}
                    c[0] = ctx.one
                    c[1] = 2*(B - A + ctx.mpq_1_4*(3*A2+B3-2)*(A2-B3) - ctx.mpq_3_16)
                    c[2] = ctx.mpq_1_2*c[1]**2 + ctx.mpq_1_16*(-16*(2*A2-3)*(B-A) + 32*R +\
                        4*(-8*A2**2 + 11*A2 + 8*A + B3 - 2)*(A2-B3)-3)
                    s1 = 0
                    s2 = 0
                    k = 0
                    tprev = 0
                    while 1:
                        if k not in c:
                            uu1 = (k-2*X-3)*(k-2*X-2*b1-1)*(k-2*X-2*b2-1)*\
                                (k-2*X-2*b3-1)
                            uu2 = (4*(k-1)**3 - 6*(4*X+B3)*(k-1)**2 + \
                                2*(24*X**2+12*B3*X+4*B+B3-1)*(k-1) - 32*X**3 - \
                                24*B3*X**2 - 4*B - 8*R - 4*(4*B+B3-1)*X + 2*B3-1)
                            uu3 = (5*(k-1)**2+2*(-10*X+A2-3*B3+3)*(k-1)+2*c[1])
                            c[k] = ctx.one/(2*k)*(uu1*c[k-3]-uu2*c[k-2]+uu3*c[k-1])
                        w = c[k] * ctx.power(-z, -0.5*k)
                        t1 = (-ctx.j)**k * ctx.mpf(2)**(-k) * w
                        t2 = ctx.j**k * ctx.mpf(2)**(-k) * w
                        if abs(t1) < 0.1*ctx.eps:
                            break
                        # Quit if the series doesn't converge quickly enough
                        if k > 5 and abs(tprev) / abs(t1) < 1.5:
                            raise ctx.NoConvergence
                        s1 += t1
                        s2 += t2
                        tprev = t1
                        k += 1
                    S = ctx.expj(ctx.pi*X+2*ctx.sqrt(-z))*s1 + \
                        ctx.expj(-(ctx.pi*X+2*ctx.sqrt(-z)))*s2
                    T1 = [0.5*S, ctx.pi, -z], [1, -0.5, X], [b1, b2, b3], [a1, a2],\
                        [], [], 0
                    T2 = [-z], [-a1], [b1,b2,b3,a2-a1],[a2,b1-a1,b2-a1,b3-a1], \
                        [a1,a1-b1+1,a1-b2+1,a1-b3+1], [a1-a2+1], 1/z
                    T3 = [-z], [-a2], [b1,b2,b3,a1-a2],[a1,b1-a2,b2-a2,b3-a2], \
                        [a2,a2-b1+1,a2-b2+1,a2-b3+1],[-a1+a2+1], 1/z
                    return T1, T2, T3
                v = ctx.hypercomb(h, [a1,a2,b1,b2,b3], force_series=True, maxterms=4*ctx.prec)
                if sum(ctx._is_real_type(u) for u in [a1,a2,b1,b2,b3,z]) == 6:
                    v = ctx.re(v)
                return v
            except ctx.NoConvergence:
                pass
        finally:
            ctx.prec = orig

    return ctx.hypsum(2, 3, (a1type, a2type, b1type, b2type, b3type), [a1, a2, b1, b2, b3], z, **kwargs)

@defun
def _hyp2f0(ctx, a_s, b_s, z, **kwargs):
    (a, atype), (b, btype) = a_s
    # We want to try aggressively to use the asymptotic expansion,
    # and fall back only when absolutely necessary
    try:
        kwargsb = kwargs.copy()
        kwargsb['maxterms'] = kwargsb.get('maxterms', ctx.prec)
        return ctx.hypsum(2, 0, (atype,btype), [a,b], z, **kwargsb)
    except ctx.NoConvergence:
        if kwargs.get('force_series'):
            raise
        pass
    def h(a, b):
        w = ctx.sinpi(b)
        rz = -1/z
        T1 = ([ctx.pi,w,rz],[1,-1,a],[],[a-b+1,b],[a],[b],rz)
        T2 = ([-ctx.pi,w,rz],[1,-1,1+a-b],[],[a,2-b],[a-b+1],[2-b],rz)
        return T1, T2
    return ctx.hypercomb(h, [a, 1+a-b], **kwargs)

@defun
def hyperu(ctx, a, b, z, **kwargs):
    a, atype = ctx._convert_param(a)
    b, btype = ctx._convert_param(b)
    z = ctx.convert(z)
    if not z:
        if ctx.re(b) <= 1:
            return ctx.gammaprod([1-b],[a-b+1])
        else:
            return ctx.inf + z
    bb = 1+a-b
    bb, bbtype = ctx._convert_param(bb)
    try:
        orig = ctx.prec
        try:
            ctx.prec += 10
            v = ctx.hypsum(2, 0, (atype, bbtype), [a, bb], -1/z, maxterms=ctx.prec)
            return v / z**a
        finally:
            ctx.prec = orig
    except ctx.NoConvergence:
        pass
    def h(a,b):
        w = ctx.sinpi(b)
        T1 = ([ctx.pi,w],[1,-1],[],[a-b+1,b],[a],[b],z)
        T2 = ([-ctx.pi,w,z],[1,-1,1-b],[],[a,2-b],[a-b+1],[2-b],z)
        return T1, T2
    return ctx.hypercomb(h, [a,b], **kwargs)

@defun_wrapped
def _erf_complex(ctx, z):
    z2 = ctx.square_exp_arg(z, -1)
    #z2 = -z**2
    v = (2/ctx.sqrt(ctx.pi))*z * ctx.hyp1f1((1,2),(3,2), z2)
    if not ctx._re(z):
        v = ctx._im(v)*ctx.j
    return v

@defun_wrapped
def _erfc_complex(ctx, z):
    if ctx.re(z) > 2:
        z2 = ctx.square_exp_arg(z)
        nz2 = ctx.fneg(z2, exact=True)
        v = ctx.exp(nz2)/ctx.sqrt(ctx.pi) * ctx.hyperu((1,2),(1,2), z2)
    else:
        v = 1 - ctx._erf_complex(z)
    if not ctx._re(z):
        v = 1+ctx._im(v)*ctx.j
    return v

@defun
def erf(ctx, z):
    z = ctx.convert(z)
    if ctx._is_real_type(z):
        try:
            return ctx._erf(z)
        except NotImplementedError:
            pass
    if ctx._is_complex_type(z) and not z.imag:
        try:
            return type(z)(ctx._erf(z.real))
        except NotImplementedError:
            pass
    return ctx._erf_complex(z)

@defun
def erfc(ctx, z):
    z = ctx.convert(z)
    if ctx._is_real_type(z):
        try:
            return ctx._erfc(z)
        except NotImplementedError:
            pass
    if ctx._is_complex_type(z) and not z.imag:
        try:
            return type(z)(ctx._erfc(z.real))
        except NotImplementedError:
            pass
    return ctx._erfc_complex(z)

@defun
def square_exp_arg(ctx, z, mult=1, reciprocal=False):
    prec = ctx.prec*4+20
    if reciprocal:
        z2 = ctx.fmul(z, z, prec=prec)
        z2 = ctx.fdiv(ctx.one, z2, prec=prec)
    else:
        z2 = ctx.fmul(z, z, prec=prec)
    if mult != 1:
        z2 = ctx.fmul(z2, mult, exact=True)
    return z2

@defun_wrapped
def erfi(ctx, z):
    if not z:
        return z
    z2 = ctx.square_exp_arg(z)
    v = (2/ctx.sqrt(ctx.pi)*z) * ctx.hyp1f1((1,2), (3,2), z2)
    if not ctx._re(z):
        v = ctx._im(v)*ctx.j
    return v

@defun_wrapped
def erfinv(ctx, x):
    xre = ctx._re(x)
    if (xre != x) or (xre < -1) or (xre > 1):
        return ctx.bad_domain("erfinv(x) is defined only for -1 <= x <= 1")
    x = xre
    #if ctx.isnan(x): return x
    if not x: return x
    if x == 1: return ctx.inf
    if x == -1: return ctx.ninf
    if abs(x) < 0.9:
        a = 0.53728*x**3 + 0.813198*x
    else:
        # An asymptotic formula
        u = ctx.ln(2/ctx.pi/(abs(x)-1)**2)
        a = ctx.sign(x) * ctx.sqrt(u - ctx.ln(u))/ctx.sqrt(2)
    ctx.prec += 10
    return ctx.findroot(lambda t: ctx.erf(t)-x, a)

@defun_wrapped
def npdf(ctx, x, mu=0, sigma=1):
    sigma = ctx.convert(sigma)
    return ctx.exp(-(x-mu)**2/(2*sigma**2)) / (sigma*ctx.sqrt(2*ctx.pi))

@defun_wrapped
def ncdf(ctx, x, mu=0, sigma=1):
    a = (x-mu)/(sigma*ctx.sqrt(2))
    if a < 0:
        return ctx.erfc(-a)/2
    else:
        return (1+ctx.erf(a))/2

@defun_wrapped
def betainc(ctx, a, b, x1=0, x2=1, regularized=False):
    if x1 == x2:
        v = 0
    elif not x1:
        if x1 == 0 and x2 == 1:
            v = ctx.beta(a, b)
        else:
            v = x2**a * ctx.hyp2f1(a, 1-b, a+1, x2) / a
    else:
        m, d = ctx.nint_distance(a)
        if m <= 0:
            if d < -ctx.prec:
                h = +ctx.eps
                ctx.prec *= 2
                a += h
            elif d < -4:
                ctx.prec -= d
        s1 = x2**a * ctx.hyp2f1(a,1-b,a+1,x2)
        s2 = x1**a * ctx.hyp2f1(a,1-b,a+1,x1)
        v = (s1 - s2) / a
    if regularized:
        v /= ctx.beta(a,b)
    return v

@defun
def gammainc(ctx, z, a=0, b=None, regularized=False):
    regularized = bool(regularized)
    z = ctx.convert(z)
    if a is None:
        a = ctx.zero
        lower_modified = False
    else:
        a = ctx.convert(a)
        lower_modified = a != ctx.zero
    if b is None:
        b = ctx.inf
        upper_modified = False
    else:
        b = ctx.convert(b)
        upper_modified = b != ctx.inf
    # Complete gamma function
    if not (upper_modified or lower_modified):
        if regularized:
            if ctx.re(z) < 0:
                return ctx.inf
            elif ctx.re(z) > 0:
                return ctx.one
            else:
                return ctx.nan
        return ctx.gamma(z)
    if a == b:
        return ctx.zero
    # Standardize
    if ctx.re(a) > ctx.re(b):
        return -ctx.gammainc(z, b, a, regularized)
    # Generalized gamma
    if upper_modified and lower_modified:
        return +ctx._gamma3(z, a, b, regularized)
    # Upper gamma
    elif lower_modified:
        return ctx._upper_gamma(z, a, regularized)
    # Lower gamma
    elif upper_modified:
        return ctx._lower_gamma(z, b, regularized)

@defun
def _lower_gamma(ctx, z, b, regularized=False):
    # Pole
    if ctx.isnpint(z):
        return type(z)(ctx.inf)
    G = [z] * regularized
    negb = ctx.fneg(b, exact=True)
    def h(z):
        T1 = [ctx.exp(negb), b, z], [1, z, -1], [], G, [1], [1+z], b
        return (T1,)
    return ctx.hypercomb(h, [z])

@defun
def _upper_gamma(ctx, z, a, regularized=False):
    # Fast integer case, when available
    if ctx.isint(z):
        try:
            if regularized:
                # Gamma pole
                if ctx.isnpint(z):
                    return type(z)(ctx.zero)
                orig = ctx.prec
                try:
                    ctx.prec += 10
                    return ctx._gamma_upper_int(z, a) / ctx.gamma(z)
                finally:
                    ctx.prec = orig
            else:
                return ctx._gamma_upper_int(z, a)
        except NotImplementedError:
            pass
    nega = ctx.fneg(a, exact=True)
    G = [z] * regularized
    # Use 2F0 series when possible; fall back to lower gamma representation
    try:
        def h(z):
            r = z-1
            return [([ctx.exp(nega), a], [1, r], [], G, [1, -r], [], 1/nega)]
        return ctx.hypercomb(h, [z], force_series=True)
    except ctx.NoConvergence:
        def h(z):
            T1 = [], [1, z-1], [z], G, [], [], 0
            T2 = [-ctx.exp(nega), a, z], [1, z, -1], [], G, [1], [1+z], a
            return T1, T2
        return ctx.hypercomb(h, [z])

@defun
def _gamma3(ctx, z, a, b, regularized=False):
    pole = ctx.isnpint(z)
    if regularized and pole:
        return ctx.zero
    try:
        ctx.prec += 15
        # We don't know in advance whether it's better to write as a difference
        # of lower or upper gamma functions, so try both
        T1 = ctx.gammainc(z, a, regularized=regularized)
        T2 = ctx.gammainc(z, b, regularized=regularized)
        R = T1 - T2
        if ctx.mag(R) - max(ctx.mag(T1), ctx.mag(T2)) > -10:
            return R
        if not pole:
            T1 = ctx.gammainc(z, 0, b, regularized=regularized)
            T2 = ctx.gammainc(z, 0, a, regularized=regularized)
            R = T1 - T2
            # May be ok, but should probably at least print a warning
            # about possible cancellation
            if 1: #ctx.mag(R) - max(ctx.mag(T1), ctx.mag(T2)) > -10:
                return R
    finally:
        ctx.prec -= 15
    raise NotImplementedError

@defun_wrapped
def expint(ctx, n, z):
    if ctx.isint(n) and ctx._is_real_type(z):
        try:
            return ctx._expint_int(n, z)
        except NotImplementedError:
            pass
    if ctx.isnan(n) or ctx.isnan(z):
        return z*n
    if z == ctx.inf:
        return 1/z
    if z == 0:
        # integral from 1 to infinity of t^n
        if ctx.re(n) <= 1:
            # TODO: reasonable sign of infinity
            return type(z)(ctx.inf)
        else:
            return ctx.one/(n-1)
    if n == 0:
        return ctx.exp(-z)/z
    if n == -1:
        return ctx.exp(-z)*(z+1)/z**2
    return z**(n-1) * ctx.gammainc(1-n, z)

@defun_wrapped
def li(ctx, z, offset=False):
    if offset:
        if z == 2:
            return ctx.zero
        return ctx.ei(ctx.ln(z)) - ctx.ei(ctx.ln2)
    if not z:
        return z
    if z == 1:
        return ctx.ninf
    return ctx.ei(ctx.ln(z))

@defun
def ei(ctx, z):
    try:
        return ctx._ei(z)
    except NotImplementedError:
        return ctx._ei_generic(z)

@defun_wrapped
def _ei_generic(ctx, z):
    # Note: the following is currently untested because mp and fp
    # both use special-case ei code
    if z == ctx.inf:
        return z
    if z == ctx.ninf:
        return ctx.zero
    if ctx.mag(z) > 1:
        try:
            r = ctx.one/z
            v = ctx.exp(z)*ctx.hyper([1,1],[],r,
                maxterms=ctx.prec, force_series=True)/z
            im = ctx._im(z)
            if im > 0:
                v += ctx.pi*ctx.j
            if im < 0:
                v -= ctx.pi*ctx.j
            return v
        except ctx.NoConvergence:
            pass
    v = z*ctx.hyp2f2(1,1,2,2,z) + ctx.euler
    if ctx._im(z):
        v += 0.5*(ctx.log(z) - ctx.log(ctx.one/z))
    else:
        v += ctx.log(abs(z))
    return v

@defun
def e1(ctx, z):
    try:
        return ctx._e1(z)
    except NotImplementedError:
        return ctx.expint(1, z)

@defun
def ci(ctx, z):
    try:
        return ctx._ci(z)
    except NotImplementedError:
        return ctx._ci_generic(z)

@defun_wrapped
def _ci_generic(ctx, z):
    if ctx.isinf(z):
        if z == ctx.inf: return ctx.zero
        if z == ctx.ninf: return ctx.pi*1j
    jz = ctx.fmul(ctx.j,z,exact=True)
    njz = ctx.fneg(jz,exact=True)
    v = 0.5*(ctx.ei(jz) + ctx.ei(njz))
    zreal = ctx._re(z)
    zimag = ctx._im(z)
    if zreal == 0:
        if zimag > 0: v += ctx.pi*0.5j
        if zimag < 0: v -= ctx.pi*0.5j
    if zreal < 0:
        if zimag >= 0: v += ctx.pi*1j
        if zimag <  0: v -= ctx.pi*1j
    if ctx._is_real_type(z) and zreal > 0:
        v = ctx._re(v)
    return v

@defun
def si(ctx, z):
    try:
        return ctx._si(z)
    except NotImplementedError:
        return ctx._si_generic(z)

@defun_wrapped
def _si_generic(ctx, z):
    if ctx.isinf(z):
        if z == ctx.inf: return 0.5*ctx.pi
        if z == ctx.ninf: return -0.5*ctx.pi
    # Suffers from cancellation near 0
    if ctx.mag(z) >= -1:
        jz = ctx.fmul(ctx.j,z,exact=True)
        njz = ctx.fneg(jz,exact=True)
        v = (-0.5j)*(ctx.ei(jz) - ctx.ei(njz))
        zreal = ctx._re(z)
        if zreal > 0:
            v -= 0.5*ctx.pi
        if zreal < 0:
            v += 0.5*ctx.pi
        if ctx._is_real_type(z):
            v = ctx._re(v)
        return v
    else:
        return z*ctx.hyp1f2((1,2),(3,2),(3,2),-0.25*z*z)

@defun_wrapped
def chi(ctx, z):
    nz = ctx.fneg(z, exact=True)
    v = 0.5*(ctx.ei(z) + ctx.ei(nz))
    zreal = ctx._re(z)
    zimag = ctx._im(z)
    if zimag > 0:
        v += ctx.pi*0.5j
    elif zimag < 0:
        v -= ctx.pi*0.5j
    elif zreal < 0:
        v += ctx.pi*1j
    return v

@defun_wrapped
def shi(ctx, z):
    # Suffers from cancellation near 0
    if ctx.mag(z) >= -1:
        nz = ctx.fneg(z, exact=True)
        v = 0.5*(ctx.ei(z) - ctx.ei(nz))
        zimag = ctx._im(z)
        if zimag > 0: v -= 0.5j*ctx.pi
        if zimag < 0: v += 0.5j*ctx.pi
        return v
    else:
        return z * ctx.hyp1f2((1,2),(3,2),(3,2),0.25*z*z)

@defun_wrapped
def fresnels(ctx, z):
    if z == ctx.inf:
        return ctx.mpf(0.5)
    if z == ctx.ninf:
        return ctx.mpf(-0.5)
    return ctx.pi*z**3/6*ctx.hyp1f2((3,4),(3,2),(7,4),-ctx.pi**2*z**4/16)

@defun_wrapped
def fresnelc(ctx, z):
    if z == ctx.inf:
        return ctx.mpf(0.5)
    if z == ctx.ninf:
        return ctx.mpf(-0.5)
    return z*ctx.hyp1f2((1,4),(1,2),(5,4),-ctx.pi**2*z**4/16)

@defun_wrapped
def airyai(ctx, z):
    if z == ctx.inf or z == ctx.ninf:
        return ctx.zero
    if z:
        # Account for exponential scaling
        ctx.prec += max(0, int(1.5*ctx.mag(z)))
    if ctx._re(z) > 4:
        # We could still use 1F1, but it results in huge cancellation;
        # the following expansion is better
        w = z**1.5
        r = -ctx.mpf(3)/(4*w)
        v = ctx.exp(-2*w/3)/(2*ctx.sqrt(ctx.pi)*ctx.nthroot(z,4))
        v *= ctx.hyp2f0((1,6),(5,6),r)
        return v
    elif ctx._re(z) > 1:
        # If not using asymptotic series:
        # cancellation: both terms are ~ 2^(z^1.5),
        # result is ~ 2^(-z^1.5), so need ~2*z^1.5 extra bits
        ctx.prec += 2*int(ctx._re(z)**1.5)
    z3 = z**3 / 9
    a = ctx.hyp0f1((2,3), z3) / (ctx.cbrt(9) * ctx.gamma(ctx.mpf(2)/3))
    b = z * ctx.hyp0f1((4,3), z3) / (ctx.cbrt(3) * ctx.gamma(ctx.mpf(1)/3))
    return a - b

@defun_wrapped
def airybi(ctx, z):
    if z == ctx.inf:
        return z
    if z == ctx.ninf:
        return 1/z
    if z:
        # Account for exponential scaling
        ctx.prec += max(0, int(1.5*ctx.mag(z)))
    z3 = z**3 / 9
    rt = ctx.nthroot(3, 6)
    a = ctx.hyp0f1((2,3), z3) / (rt * ctx.gamma(ctx.mpf(2)/3))
    b = z * rt * ctx.hyp0f1((4,3), z3) / ctx.gamma(ctx.mpf(1)/3)
    return a + b

@defun_wrapped
def hermite(ctx, n, z, **kwargs):
    if not z:
        try:
            return 2**n * ctx.sqrt(ctx.pi) / ctx.gamma(0.5*(1-n))
        except ValueError:
            return 0.0*(n+z)
    if ctx.re(z) > 0 or (ctx.re(z) == 0 and ctx.im(z) > 0) or ctx.isnpint(-n):
        prec = ctx.prec
        ctx.prec = ctx.prec*4+20
        z2 = -z**(-2)
        ctx.prec = prec
        return (2*z)**n * ctx.hyp2f0(-0.5*n, -0.5*(n-1), z2, **kwargs)
    else:
        prec = ctx.prec
        ctx.prec = ctx.prec*4+20
        z2 = z**2
        ctx.prec = prec
        return ctx.hermite(n,-z) + 2**(n+2)*ctx.sqrt(ctx.pi) * (-z) / \
            ctx.gamma(-0.5*n) * ctx.hyp1f1((1-n)*0.5, 1.5, z2, **kwargs)

@defun_wrapped
def gegenbauer(ctx, n, a, z, **kwargs):
    # Special cases: a+0.5, a*2 poles
    if ctx.isnpint(a):
        return 0*(z+n)
    if ctx.isnpint(a+0.5):
        # TODO: something else is required here
        # E.g.: gegenbauer(-2, -0.5, 3) == -12
        if ctx.isnpint(n+1):
            raise NotImplementedError("Gegenbauer function with two limits")
        def h(a):
            a2 = 2*a
            T = [], [], [n+a2], [n+1, a2], [-n, n+a2], [a+0.5], 0.5*(1-z)
            return [T]
        return ctx.hypercomb(h, [a], **kwargs)
    def h(n):
        a2 = 2*a
        T = [], [], [n+a2], [n+1, a2], [-n, n+a2], [a+0.5], 0.5*(1-z)
        return [T]
    return ctx.hypercomb(h, [n], **kwargs)

@defun_wrapped
def jacobi(ctx, n, a, b, x, **kwargs):
    if not ctx.isnpint(a):
        def h(n):
            return (([], [], [a+n+1], [n+1, a+1], [-n, a+b+n+1], [a+1], (1-x)*0.5),)
        return ctx.hypercomb(h, [n], **kwargs)
    if not ctx.isint(b):
        def h(n, a):
            return (([], [], [-b], [n+1, -b-n], [-n, a+b+n+1], [b+1], (x+1)*0.5),)
        return ctx.hypercomb(h, [n, a], **kwargs)
    # XXX: determine appropriate limit
    return ctx.binomial(n+a,n) * ctx.hyp2f1(-n,1+n+a+b,a+1,(1-x)/2, **kwargs)

@defun_wrapped
def laguerre(ctx, n, a, z, **kwargs):
    # XXX: limits, poles
    #if ctx.isnpint(n):
    #    return 0*(a+z)
    def h(a):
        return (([], [], [a+n+1], [a+1, n+1], [-n], [a+1], z),)
    return ctx.hypercomb(h, [a], **kwargs)

@defun_wrapped
def legendre(ctx, n, x, **kwargs):
    if ctx.isint(n):
        n = int(n)
        # Accuracy near zeros
        if (n + (n < 0)) & 1:
            if not x:
                return x
            mag = ctx.mag(x)
            if mag < -2*ctx.prec-10:
                return x
            if mag < -5:
                ctx.prec += -mag
    return ctx.hyp2f1(-n,n+1,1,(1-x)/2, **kwargs)

@defun
def legenp(ctx, n, m, z, type=2, **kwargs):
    # Legendre function, 1st kind
    n = ctx.convert(n)
    m = ctx.convert(m)
    # Faster
    if not m:
        return ctx.legendre(n, z, **kwargs)
    # TODO: correct evaluation at singularities
    if type == 2:
        def h(n,m):
            g = m*0.5
            T = [1+z, 1-z], [g, -g], [], [1-m], [-n, n+1], [1-m], 0.5*(1-z)
            return (T,)
        return ctx.hypercomb(h, [n,m], **kwargs)
    if type == 3:
        def h(n,m):
            g = m*0.5
            T = [z+1, z-1], [g, -g], [], [1-m], [-n, n+1], [1-m], 0.5*(1-z)
            return (T,)
        return ctx.hypercomb(h, [n,m], **kwargs)
    raise ValueError("requires type=2 or type=3")

@defun
def legenq(ctx, n, m, z, type=2, **kwargs):
    # Legendre function, 2nd kind
    n = ctx.convert(n)
    m = ctx.convert(m)
    z = ctx.convert(z)
    if z in (1, -1):
        #if ctx.isint(m):
        #    return ctx.nan
        #return ctx.inf  # unsigned
        return ctx.nan
    if type == 2:
        def h(n, m):
            s = 2 * ctx.sinpi(m) / ctx.pi
            c = ctx.cospi(m)
            a = 1+z
            b = 1-z
            u = m/2
            w = (1-z)/2
            T1 = [s, c, a, b], [-1, 1, u, -u], [], [1-m], \
                [-n, n+1], [1-m], w
            T2 = [-s, a, b], [-1, -u, u], [n+m+1], [n-m+1, m+1], \
                [-n, n+1], [m+1], w
            return T1, T2
        return ctx.hypercomb(h, [n, m], **kwargs)
    if type == 3:
        # The following is faster when there only is a single series
        # Note: not valid for -1 < z < 0 (?)
        if abs(z) > 1:
            def h(n, m):
                T1 = [ctx.expjpi(m), 2, ctx.pi, z, z-1, z+1], \
                     [1, -n-1, 0.5, -n-m-1, 0.5*m, 0.5*m], \
                     [n+m+1], [n+1.5], \
                     [0.5*(2+n+m), 0.5*(1+n+m)], [n+1.5], z**(-2)
                return [T1]
            return ctx.hypercomb(h, [n, m], **kwargs)
        else:
            # not valid for 1 < z < inf ?
            def h(n, m):
                s = 2 * ctx.sinpi(m) / ctx.pi
                c = ctx.expjpi(m)
                a = 1+z
                b = z-1
                u = m/2
                w = (1-z)/2
                T1 = [s, c, a, b], [-1, 1, u, -u], [], [1-m], \
                    [-n, n+1], [1-m], w
                T2 = [-s, c, a, b], [-1, 1, -u, u], [n+m+1], [n-m+1, m+1], \
                    [-n, n+1], [m+1], w
                return T1, T2
            return ctx.hypercomb(h, [n, m], **kwargs)
    raise ValueError("requires type=2 or type=3")

@defun_wrapped
def chebyt(ctx, n, x, **kwargs):
    return ctx.hyp2f1(-n,n,(1,2),(1-x)/2, **kwargs)

@defun_wrapped
def chebyu(ctx, n, x, **kwargs):
    return (n+1) * ctx.hyp2f1(-n, n+2, (3,2), (1-x)/2, **kwargs)

@defun
def j0(ctx, x):
    """Computes the Bessel function `J_0(x)`. See :func:`besselj`."""
    return ctx.besselj(0, x)

@defun
def j1(ctx, x):
    """Computes the Bessel function `J_1(x)`.  See :func:`besselj`."""
    return ctx.besselj(1, x)

@defun
def besselj(ctx, n, z, derivative=0, **kwargs):
    if type(n) is int:
        n_isint = True
    else:
        n = ctx.convert(n)
        n_isint = ctx.isint(n)
        if n_isint:
            n = int(n)
    if n_isint and n < 0:
        return (-1)**n * ctx.besselj(-n, z, derivative, **kwargs)
    z = ctx.convert(z)
    M = ctx.mag(z)
    if derivative:
        d = ctx.convert(derivative)
        # TODO: the integer special-casing shouldn't be necessary.
        # However, the hypergeometric series gets inaccurate for large d
        # because of inaccurate pole cancellation at a pole far from
        # zero (needs to be fixed in hypercomb or hypsum)
        if ctx.isint(d) and d >= 0:
            d = int(d)
            orig = ctx.prec
            try:
                ctx.prec += 15
                v = ctx.fsum((-1)**k * ctx.binomial(d,k) * ctx.besselj(2*k+n-d,z)
                    for k in range(d+1))
            finally:
                ctx.prec = orig
            v *= ctx.mpf(2)**(-d)
        else:
            def h(n,d):
                r = ctx.fmul(ctx.fmul(z, z, prec=ctx.prec+M), -0.25, exact=True)
                B = [0.5*(n-d+1), 0.5*(n-d+2)]
                T = [([2,ctx.pi,z],[d-2*n,0.5,n-d],[],B,[(n+1)*0.5,(n+2)*0.5],B+[n+1],r)]
                return T
            v = ctx.hypercomb(h, [n,d], **kwargs)
    else:
    # Fast case: J_n(x), n int, appropriate magnitude for fixed-point calculation
        if (not derivative) and n_isint and abs(M) < 10 and abs(n) < 20:
            try:
                return ctx._besselj(n, z)
            except NotImplementedError:
                pass
        if not z:
            if not n:
                v = ctx.one + n+z
            elif ctx.re(n) > 0:
                v = n*z
            else:
                v = ctx.inf + z + n
        else:
            #v = 0
            orig = ctx.prec
            try:
                # XXX: workaround for accuracy in low level hypergeometric series
                # when alternating, large arguments
                ctx.prec += min(3*abs(M), ctx.prec)
                w = ctx.fmul(z, 0.5, exact=True)
                def h(n):
                    r = ctx.fneg(ctx.fmul(w, w, prec=max(0,ctx.prec+M)), exact=True)
                    return [([w], [n], [], [n+1], [], [n+1], r)]
                v = ctx.hypercomb(h, [n], **kwargs)
            finally:
                ctx.prec = orig
        v = +v
    return v

@defun
def besseli(ctx, n, z, derivative=0, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    if not z:
        if derivative:
            raise ValueError
        if not n:
            # I(0,0) = 1
            return 1+n+z
        if ctx.isint(n):
            return 0*(n+z)
        r = ctx.re(n)
        if r == 0:
            return ctx.nan*(n+z)
        elif r > 0:
            return 0*(n+z)
        else:
            return ctx.inf+(n+z)
    M = ctx.mag(z)
    if derivative:
        d = ctx.convert(derivative)
        def h(n,d):
            r = ctx.fmul(ctx.fmul(z, z, prec=ctx.prec+M), 0.25, exact=True)
            B = [0.5*(n-d+1), 0.5*(n-d+2), n+1]
            T = [([2,ctx.pi,z],[d-2*n,0.5,n-d],[n+1],B,[(n+1)*0.5,(n+2)*0.5],B,r)]
            return T
        v = ctx.hypercomb(h, [n,d], **kwargs)
    else:
        def h(n):
            w = ctx.fmul(z, 0.5, exact=True)
            r = ctx.fmul(w, w, prec=max(0,ctx.prec+M))
            return [([w], [n], [], [n+1], [], [n+1], r)]
        v = ctx.hypercomb(h, [n], **kwargs)
    return v

@defun_wrapped
def bessely(ctx, n, z, derivative=0, **kwargs):
    if not z:
        if derivative:
            # Not implemented
            raise ValueError
        if not n:
            # ~ log(z/2)
            return -ctx.inf + (n+z)
        if ctx.im(n):
            return nan * (n+z)
        r = ctx.re(n)
        q = n+0.5
        if ctx.isint(q):
            if n > 0:
                return -ctx.inf + (n+z)
            else:
                return 0 * (n+z)
        if r < 0 and int(ctx.floor(q)) % 2:
            return ctx.inf + (n+z)
        else:
            return ctx.ninf + (n+z)
    # XXX: use hypercomb
    ctx.prec += 10
    m, d = ctx.nint_distance(n)
    if d < -ctx.prec:
        h = +ctx.eps
        ctx.prec *= 2
        n += h
    elif d < 0:
        ctx.prec -= d
    # TODO: avoid cancellation for imaginary arguments
    return (ctx.besselj(n,z,derivative,**kwargs)*ctx.cospi(n) - \
        ctx.besselj(-n,z,derivative,**kwargs))/ctx.sinpi(n)

@defun_wrapped
def besselk(ctx, n, z, **kwargs):
    if not z:
        return ctx.inf
    M = ctx.mag(z)
    if M < 1:
        # Represent as limit definition
        def h(n):
            r = (z/2)**2
            T1 = [z, 2], [-n, n-1], [n], [], [], [1-n], r
            T2 = [z, 2], [n, -n-1], [-n], [], [], [1+n], r
            return T1, T2
    # We could use the limit definition always, but it leads
    # to very bad cancellation (of exponentially large terms)
    # for large real z
    # Instead represent in terms of 2F0
    else:
        ctx.prec += M
        def h(n):
            return [([ctx.pi/2, z, ctx.exp(-z)], [0.5,-0.5,1], [], [], \
                [n+0.5, 0.5-n], [], -1/(2*z))]
    return ctx.hypercomb(h, [n], **kwargs)

@defun_wrapped
def hankel1(ctx,n,x,**kwargs):
    return ctx.besselj(n,x,**kwargs) + ctx.j*ctx.bessely(n,x,**kwargs)

@defun_wrapped
def hankel2(ctx,n,x,**kwargs):
    return ctx.besselj(n,x,**kwargs) - ctx.j*ctx.bessely(n,x,**kwargs)

@defun_wrapped
def whitm(ctx,k,m,z,**kwargs):
    if z == 0:
        # M(k,m,z) = 0^(1/2+m)
        if ctx.re(m) > -0.5:
            return z
        elif ctx.re(m) < -0.5:
            return ctx.inf + z
        else:
            return ctx.nan * z
    x = ctx.fmul(-0.5, z, exact=True)
    y = 0.5+m
    return ctx.exp(x) * z**y * ctx.hyp1f1(y-k, 1+2*m, z, **kwargs)

@defun_wrapped
def whitw(ctx,k,m,z,**kwargs):
    if z == 0:
        g = abs(ctx.re(m))
        if g < 0.5:
            return z
        elif g > 0.5:
            return ctx.inf + z
        else:
            return ctx.nan * z
    x = ctx.fmul(-0.5, z, exact=True)
    y = 0.5+m
    return ctx.exp(x) * z**y * ctx.hyperu(y-k, 1+2*m, z, **kwargs)

@defun
def struveh(ctx,n,z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/StruveH/26/01/02/
    def h(n):
        return [([z/2, 0.5*ctx.sqrt(ctx.pi)], [n+1, -1], [], [n+1.5], [1], [1.5, n+1.5], -(z/2)**2)]
    return ctx.hypercomb(h, [n], **kwargs)

@defun
def struvel(ctx,n,z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/StruveL/26/01/02/
    def h(n):
        return [([z/2, 0.5*ctx.sqrt(ctx.pi)], [n+1, -1], [], [n+1.5], [1], [1.5, n+1.5], (z/2)**2)]
    return ctx.hypercomb(h, [n], **kwargs)

@defun
def ber(ctx, n, z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinBer2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        T1 = [ctx.cospi(0.75*n), z/2], [1, n], [], [n+1], [], [0.5, 0.5*(n+1), 0.5*n+1], r
        T2 = [-ctx.sinpi(0.75*n), z/2], [1, n+2], [], [n+2], [], [1.5, 0.5*(n+3), 0.5*n+1], r
        return T1, T2
    return ctx.hypercomb(h, [n], **kwargs)

@defun
def bei(ctx, n, z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinBei2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        T1 = [ctx.cospi(0.75*n), z/2], [1, n+2], [], [n+2], [], [1.5, 0.5*(n+3), 0.5*n+1], r
        T2 = [ctx.sinpi(0.75*n), z/2], [1, n], [], [n+1], [], [0.5, 0.5*(n+1), 0.5*n+1], r
        return T1, T2
    return ctx.hypercomb(h, [n], **kwargs)

@defun
def ker(ctx, n, z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinKer2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        T1 = [2, z, 4*ctx.cospi(0.25*n)], [-n-3, n, 1], [-n], [], [], [0.5, 0.5*(1+n), 0.5*(n+2)], r
        T2 = [2, z, -ctx.sinpi(0.25*n)], [-n-3, 2+n, 1], [-n-1], [], [], [1.5, 0.5*(3+n), 0.5*(n+2)], r
        T3 = [2, z, 4*ctx.cospi(0.75*n)], [n-3, -n, 1], [n], [], [], [0.5, 0.5*(1-n), 1-0.5*n], r
        T4 = [2, z, -ctx.sinpi(0.75*n)], [n-3, 2-n, 1], [n-1], [], [], [1.5, 0.5*(3-n), 1-0.5*n], r
        return T1, T2, T3, T4
    return ctx.hypercomb(h, [n], **kwargs)

@defun
def kei(ctx, n, z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinKei2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        T1 = [-ctx.cospi(0.75*n), 2, z], [1, n-3, 2-n], [n-1], [], [], [1.5, 0.5*(3-n), 1-0.5*n], r
        T2 = [-ctx.sinpi(0.75*n), 2, z], [1, n-1, -n], [n], [], [], [0.5, 0.5*(1-n), 1-0.5*n], r
        T3 = [-ctx.sinpi(0.25*n), 2, z], [1, -n-1, n], [-n], [], [], [0.5, 0.5*(n+1), 0.5*(n+2)], r
        T4 = [-ctx.cospi(0.25*n), 2, z], [1, -n-3, n+2], [-n-1], [], [], [1.5, 0.5*(n+3), 0.5*(n+2)], r
        return T1, T2, T3, T4
    return ctx.hypercomb(h, [n], **kwargs)

@defun
def meijerg(ctx, a_s, b_s, z, r=1, series=None, **kwargs):
    an, ap = a_s
    bm, bq = b_s
    n = len(an)
    p = n + len(ap)
    m = len(bm)
    q = m + len(bq)
    a = an+ap
    b = bm+bq
    a = map(ctx.convert, a)
    b = map(ctx.convert, b)
    z = ctx.convert(z)
    if series is None:
        if p < q: series = 1
        if p > q: series = 2
        if p == q:
            if m+n == p and abs(z) > 1:
                series = 2
            else:
                series = 1
    if kwargs.get('verbose'):
        print "Meijer G m,n,p,q,series =", m,n,p,q,series
    if series == 1:
        def h(*args):
            a = args[:p]
            b = args[p:]
            terms = []
            for k in range(m):
                bases = [z]
                expts = [b[k]/r]
                gn = [b[j]-b[k] for j in range(m) if j != k]
                gn += [1-a[j]+b[k] for j in range(n)]
                gd = [a[j]-b[k] for j in range(n,p)]
                gd += [1-b[j]+b[k] for j in range(m,q)]
                hn = [1-a[j]+b[k] for j in range(p)]
                hd = [1-b[j]+b[k] for j in range(q) if j != k]
                hz = (-ctx.one)**(p-m-n) * z**(ctx.one/r)
                terms.append((bases, expts, gn, gd, hn, hd, hz))
            return terms
    else:
        def h(*args):
            a = args[:p]
            b = args[p:]
            terms = []
            for k in range(n):
                bases = [z]
                if r == 1:
                    expts = [a[k]-1]
                else:
                    expts = [(a[k]-1)/ctx.convert(r)]
                gn = [a[k]-a[j] for j in range(n) if j != k]
                gn += [1-a[k]+b[j] for j in range(m)]
                gd = [a[k]-b[j] for j in range(m,q)]
                gd += [1-a[k]+a[j] for j in range(n,p)]
                hn = [1-a[k]+b[j] for j in range(q)]
                hd = [1+a[j]-a[k] for j in range(p) if j != k]
                hz = (-ctx.one)**(q-m-n) / z**(ctx.one/r)
                terms.append((bases, expts, gn, gd, hn, hd, hz))
            return terms
    return ctx.hypercomb(h, a+b, **kwargs)

@defun_wrapped
def appellf1(ctx,a,b1,b2,c,z1,z2,**kwargs):
    # Assume z1 smaller
    # We will use z1 for the outer loop
    if abs(z1) > abs(z2):
        z1, z2 = z2, z1
        b1, b2 = b2, b1
    def ok(x):
        return abs(x) < 0.99
    # Finite cases
    if ctx.isnpint(a):
        pass
    elif ctx.isnpint(b1):
        pass
    elif ctx.isnpint(b2):
        z1, z2, b1, b2 = z2, z1, b2, b1
    else:
        #print z1, z2
        # Note: ok if |z2| > 1, because
        # 2F1 implements analytic continuation
        if not ok(z1):
            u1 = (z1-z2)/(z1-1)
            if not ok(u1):
                raise ValueError("Analytic continuation not implemented")
            #print "Using analytic continuation"
            return (1-z1)**(-b1)*(1-z2)**(c-a-b2)*\
                ctx.appellf1(c-a,b1,c-b1-b2,c,u1,z2,**kwargs)
    #print "inner is", a, b2, c
    one = ctx.one
    s = 0
    t = 1
    k = 0
    while 1:
        h = ctx.hyp2f1(a,b2,c,z2,zeroprec=ctx.prec,**kwargs)
        term = t * h
        if abs(term) < ctx.eps and abs(h) > 10*ctx.eps:
            break
        s += term
        k += 1
        t = (t*a*b1*z1) / (c*k)
        c += one
        a += one
        b1 += one
    return s

@defun_wrapped
def coulombc(ctx, l, eta, _cache={}):
    if (l, eta) in _cache and _cache[l,eta][0] >= ctx.prec:
        return +_cache[l,eta][1]
    G3 = ctx.loggamma(2*l+2)
    G1 = ctx.loggamma(1+l+ctx.j*eta)
    G2 = ctx.loggamma(1+l-ctx.j*eta)
    v = 2**l * ctx.exp((-ctx.pi*eta+G1+G2)/2 - G3)
    if not (ctx.im(l) or ctx.im(eta)):
        v = ctx.re(v)
    _cache[l,eta] = (ctx.prec, v)
    return v

@defun_wrapped
def coulombf(ctx, l, eta, z, w=1, chop=True, **kwargs):
    # Regular Coulomb wave function
    # Note: w can be either 1 or -1; the other may be better in some cases
    # TODO: check that chop=True chops when and only when it should
    #ctx.prec += 10
    def h(l, eta):
        try:
            jw = ctx.j*w
            jwz = ctx.fmul(jw, z, exact=True)
            jwz2 = ctx.fmul(jwz, -2, exact=True)
            C = ctx.coulombc(l, eta)
            T1 = [C, z, ctx.exp(jwz)], [1, l+1, 1], [], [], [1+l+jw*eta], \
                [2*l+2], jwz2
        except ValueError:
            T1 = [0], [-1], [], [], [], [], 0
        return (T1,)
    v = ctx.hypercomb(h, [l,eta], **kwargs)
    if chop and (not ctx.im(l)) and (not ctx.im(eta)) and (not ctx.im(z)) and \
        (ctx.re(z) >= 0):
        v = ctx.re(v)
    return v

@defun_wrapped
def _coulomb_chi(ctx, l, eta, _cache={}):
    if (l, eta) in _cache and _cache[l,eta][0] >= ctx.prec:
        return _cache[l,eta][1]
    def terms():
        l2 = -l-1
        jeta = ctx.j*eta
        return [ctx.loggamma(1+l+jeta) * (-0.5j),
            ctx.loggamma(1+l-jeta) * (0.5j),
            ctx.loggamma(1+l2+jeta) * (0.5j),
            ctx.loggamma(1+l2-jeta) * (-0.5j),
            -(l+0.5)*ctx.pi]
    v = ctx.sum_accurately(terms, 1)
    _cache[l,eta] = (ctx.prec, v)
    return v

@defun_wrapped
def coulombg(ctx, l, eta, z, w=1, chop=True, **kwargs):
    # Irregular Coulomb wave function
    # Note: w can be either 1 or -1; the other may be better in some cases
    # TODO: check that chop=True chops when and only when it should
    if not ctx._im(l):
        l = ctx._re(l)  # XXX: for isint
    def h(l, eta):
        # Force perturbation for integers and half-integers
        if ctx.isint(l*2):
            T1 = [0], [-1], [], [], [], [], 0
            return (T1,)
        l2 = -l-1
        try:
            chi = ctx._coulomb_chi(l, eta)
            jw = ctx.j*w
            s = ctx.sin(chi); c = ctx.cos(chi)
            C1 = ctx.coulombc(l,eta)
            C2 = ctx.coulombc(l2,eta)
            u = ctx.exp(jw*z)
            x = -2*jw*z
            T1 = [s, C1, z, u, c], [-1, 1, l+1, 1, 1], [], [], \
                [1+l+jw*eta], [2*l+2], x
            T2 = [-s, C2, z, u],   [-1, 1, l2+1, 1],    [], [], \
                [1+l2+jw*eta], [2*l2+2], x
            return T1, T2
        except ValueError:
            T1 = [0], [-1], [], [], [], [], 0
            return (T1,)
    v = ctx.hypercomb(h, [l,eta], **kwargs)
    if chop and (not ctx._im(l)) and (not ctx._im(eta)) and (not ctx._im(z)) and \
        (ctx._re(z) >= 0):
        v = ctx._re(v)
    return v

@defun
def spherharm(ctx, l, m, theta, phi, **kwargs):
    l = ctx.convert(l)
    m = ctx.convert(m)
    theta = ctx.convert(theta)
    phi = ctx.convert(phi)
    l_isint = ctx.isint(l)
    l_natural = l_isint and l >= 0
    m_isint = ctx.isint(m)
    if l_isint and l < 0 and m_isint:
        return ctx.spherharm(-(l+1), m, theta, phi, **kwargs)
    if theta == 0 and m_isint and m < 0:
        return ctx.zero * 1j
    if l_natural and m_isint:
        if abs(m) > l:
            return ctx.zero * 1j
        # http://functions.wolfram.com/Polynomials/
        #     SphericalHarmonicY/26/01/02/0004/
        def h(l,m):
            C = [-1, ctx.expj(m*phi),
                 (2*l+1)*ctx.fac(l+abs(m))/ctx.pi/ctx.fac(l-abs(m)),
                 ctx.sin(theta)**2,
                 ctx.fac(abs(m)), 2]
            P = [0.5*m*(ctx.sign(m)+1), 1, 0.5, 0.5*abs(m), -1, -abs(m)-1]
            return ((C, P, [], [], [abs(m)-l, l+abs(m)+1], [abs(m)+1],
                ctx.sin(0.5*theta)**2),)
    else:
        # http://functions.wolfram.com/HypergeometricFunctions/
        #     SphericalHarmonicYGeneral/26/01/02/0001/
        def h(l,m):
            if ctx.isnpint(l-m+1) or ctx.isnpint(l+m+1) or ctx.isnpint(1-m):
                return (([0], [-1], [], [], [], [], 0),)
            C = [0.5*ctx.expj(m*phi),
                 (2*l+1)/ctx.pi,
                 ctx.gamma(l-m+1),
                 ctx.gamma(l+m+1),
                 ctx.cos(0.5*theta)**2,
                 ctx.sin(0.5*theta)**2]
            P = [1, 0.5, 0.5, -0.5, 0.5*m, -0.5*m]
            return ((C, P, [], [1-m], [-l,l+1], [1-m], ctx.sin(0.5*theta)**2),)
    return ctx.hypercomb(h, [l,m], **kwargs)
