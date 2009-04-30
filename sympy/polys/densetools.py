"""Advanced tools for dense recursive polynomials in `K[x]` or `K[X]`. """

from sympy.polys.densebasic import (
    dup_strip, dmp_strip,
    dup_convert, dmp_convert,
    dup_degree, dmp_degree,
    dup_LC, dmp_LC, dmp_ground_LC,
    dup_TC, dmp_TC, dmp_ground_TC,
    dmp_zero, dmp_one, dmp_ground,
    dmp_zero_p, dmp_one_p,
    dmp_multi_deflate, dmp_inflate,
    dup_to_raw_dict, dup_from_raw_dict,
    dmp_zeros,
)

from sympy.polys.densearith import (
    dup_add_term,
    dmp_mul_term,
    dup_neg, dmp_neg,
    dup_add, dmp_add,
    dup_sub, dmp_sub,
    dup_mul, dmp_mul,
    dup_pow, dmp_pow,
    dup_div, dmp_div,
    dup_rem, dmp_rem,
    dup_quo, dmp_quo,
    dup_exquo, dmp_exquo,
    dup_prem, dmp_prem,
    dup_add_mul, dup_sub_mul,
    dup_mul_ground, dmp_mul_ground,
    dup_quo_ground, dmp_quo_ground,
    dup_exquo_ground, dmp_exquo_ground,
    dup_max_norm, dmp_max_norm,
)

from sympy.polys.polyerrors import (
    HeuristicGCDFailed,
    NotInvertible,
    DomainError,
)

def dup_ground_to_ring(f, K0, K1):
    """Clear denominators, i.e. transform `K_0` to `K_1`, but don't normalize. """
    common = K1.one

    for c in f:
        common = K1.lcm(common, K0.denom(c))

    if K1.is_one(common):
        return common, f
    else:
        return common, dup_mul_ground(f, common, K0)

def dmp_ground_to_ring(f, u, K0, K1):
    """Clear denominators, i.e. transform `K_0` to `K_1`, but don't normalize. """
    if not u:
        return dup_ground_to_ring(f, K0, K1)

    def rec_ground_to_ring(g, v):
        common = K1.one

        if not v:
            for c in g:
                common = K1.lcm(common, K0.denom(c))
        else:
            for c in g:
                common = K1.lcm(common, rec_ground_to_ring(c, v-1))

        return common

    common = rec_ground_to_ring(f, u)

    if not K1.is_one(common):
        f = dmp_mul_ground(f, common, u, K0)

    return common, f

def dup_integrate(f, m, K):
    """Computes indefinite integral of `f` in `K[x]`. """
    if m <= 0 or not f:
        return f

    g = [K.zero]*m

    for i, c in enumerate(reversed(f)):
        n = i+1

        for j in xrange(1, m):
            n *= i+j+1

        g.insert(0, K.quo(c, n))

    return g

def dmp_integrate(f, m, u, K):
    """Computes indefinite integral of `f` in `x_0` in `K[X]`. """
    if not u:
        return dup_integrate(f, m, K)

    if m <= 0 or dmp_zero_p(f, u):
        return f

    g = dmp_zeros(m, u-1, K)

    for i, c in enumerate(reversed(f)):
        n = i+1

        for j in xrange(1, m):
            n *= i+j+1

        g.insert(0, dmp_quo_ground(c, K(n), u-1, K))

    return g

def dmp_integrate_in(f, m, j, u, K):
    """Computes indefinite integral of `f` in `x_j` in `K[X]`. """
    if j < 0 or j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))

    def rec_integrate_in(g, v, i):
        if i == j:
            return dmp_integrate(g, m, v, K)
        else:
            return dmp_strip([ rec_integrate_in(c, v-1, i+1) for c in g ], v)

    return rec_integrate_in(f, u, 0)

def dup_diff(f, m, K):
    """m-th order derivative of a polynomial in `K[x]`. """
    if m <= 0:
        return f

    n = dup_degree(f)

    if n < m:
        return []

    deriv, c = [], K.one

    for i in xrange(0, m):
        c, n = c*n, n-1

    for coeff in f[:-m]:
        deriv.append(coeff*c)
        c, n = n*K.exquo(c, n+m), n-1

    return deriv

def dmp_diff(f, m, u, K):
    """m-th order derivative in `x_0` of a polynomial in `K[X]`. """
    if not u:
        return dup_diff(f, m, K)
    if m <= 0:
        return f

    n = dmp_degree(f, u)

    if n < m:
        return dmp_zero(u)

    deriv, c, v = [], K.one, u-1

    for i in xrange(0, m):
        c, n = c*n, n-1

    for coeff in f[:-m]:
        h = dmp_mul_ground(coeff, c, v, K)
        c, n = n*K.exquo(c, n+m), n-1
        deriv.append(h)

    return deriv

def dmp_diff_in(f, m, j, u, K):
    """m-th order derivative in `x_j` of a polynomial in `K[X]`. """
    if j < 0 or j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))

    def rec_diff_in(g, v, i):
        if i == j:
            return dmp_diff(g, m, v, K)
        else:
            return dmp_strip([ rec_diff_in(c, v-1, i+1) for c in g ], v)

    return rec_diff_in(f, u, 0)

def dup_eval(f, a, K):
    """Evaluate a polynomial at `x = a` in `K[x]` using Horner scheme. """
    if not a:
        return dup_TC(f, K)

    result = K.zero

    for c in f:
        result *= a
        result += c

    return result

def dmp_eval(f, a, u, K):
    """Evaluate a polynomial at `x_0 = a` in `K[X]` using Horner scheme. """
    if not u:
        return dup_eval(f, a, K)

    if not a:
        return dmp_TC(f, K)

    result, v = dmp_LC(f, K), u-1

    for coeff in f[1:]:
        result = dmp_mul_ground(result, a, v, K)
        result = dmp_add(result, coeff, v, K)

    return result

def dmp_eval_in(f, a, j, u, K):
    """Evaluate a polynomial at `x_j = a` in `K[X]` using Horner scheme. """
    if j < 0 or j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))

    def rec_eval_in(g, v, i):
        if i == j:
            return dmp_eval(g, a, v, K)
        else:
            return dmp_strip([ rec_eval_in(c, v-1, i+1) for c in g ], v-1)

    return rec_eval_in(f, u, 0)

def dmp_eval_tail(f, A, u, K):
    """Evaluate a polynomial at `x_j = a_j, ...` in `K[X]`. """
    def rec_eval_tail(g, i):
        if i == u:
            return dup_eval(g, A[-1], K)
        else:
            h = [ rec_eval_tail(c, i+1) for c in g ]

            if i < u - len(A) + 1:
                return h
            else:
                return dup_eval(h, A[-u+i-1], K)

    if not A:
        return f

    if dmp_zero_p(f, u):
        return dmp_zero(u - len(A))

    e = rec_eval_tail(f, 0)

    if u == len(A)-1:
        return e
    else:
        return dmp_strip(e, u - len(A))

def dmp_diff_eval_in(f, m, a, j, u, K):
    """Differentiate and evaluate a polynomial in `x_j` at `a` in `K[X]`. """
    if j < 0:
        j += u + 1

    if j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))
    if not j:
        return dmp_eval(dmp_diff(f, m, u, K), a, u, K)

    def rec_diff_eval(g, v, i):
        if i == j:
            return dmp_eval(dmp_diff(g, m, v, K), a, v, K)
        else:
            return dmp_strip([ rec_diff_eval(c, v-1, i+1) for c in g ], v-1)

    return rec_diff_eval(f, u, 0)

def dup_half_gcdex(f, g, K):
    """Half extended Euclidean algorithm in `F[x]`. """
    if not K.has_Field:
        raise DomainError('computation can be done only in a field')

    a, b = [K.one], []

    while g:
        q, r = dup_div(f, g, K)
        f, g = g, r
        a, b = b, dup_sub_mul(a, q, b, K)

    a = dup_quo_ground(a, dup_LC(f, K), K)
    f = dup_monic(f, K)

    return a, f

def dup_gcdex(f, g, K):
    """Extended Euclidean algorithm in `F[x]`. """
    s, h = dup_half_gcdex(f, g, K)

    F = dup_sub_mul(h, s, f, K)
    t = dup_exquo(F, g, K)

    return s, t, h

def dup_invert(f, g, K):
    """Compute multiplicative inverse of `f` in `F[x]/(g(x))`. """
    s, h = dup_half_gcdex(f, g, K)

    if h == [K.one]:
        return dup_rem(s, g, K)
    else:
        raise NotInvertible("zero divisor")

def dup_inner_subresultants(f, g, K):
    """Subresultant PRS algorithm in `K[x]`. """
    n = dup_degree(f)
    m = dup_degree(g)

    if n < m:
        f, g = g, f
        n, m = m, n

    R = [f, g]
    d = n - m

    b = (-K.one)**(d+1)
    c =  -K.one

    B, D = [b], [d]

    if not f or not g:
        return R, B, D

    h = dup_prem(f, g, K)
    h = dup_mul_ground(h, b, K)

    while h:
        k = dup_degree(h)
        R.append(h)

        lc = dup_LC(g, K)

        if not d:
            q = c
        else:
            q = c**(d-1)

        c = K.exquo((-lc)**d, q)
        b = -lc * c**(m-k)

        f, g, m, d = g, h, k, m-k

        B.append(b)
        D.append(d)

        h = dup_prem(f, g, K)
        h = dup_exquo_ground(h, b, K)

    return R, B, D

def dup_subresultants(f, g, K):
    """Computes subresultant PRS of two polynomials in `K[x]`. """
    return dup_inner_subresultants(f, g, K)[0]

def dup_inner_resultant(f, g, K):
    """Resultant algorithm in `K[x]` using subresultant PRS. """
    if not f or not g:
        return (K.zero, [])

    R, B, D = dup_inner_subresultants(f, g, K)

    if dup_degree(R[-1]) > 0:
        return (K.zero, R)
    if R[-2] == [K.one]:
        return (dup_LC(R[-1], K), R)

    s, i = 1, 1
    p, q = K.one, K.one

    for b, d in zip(B, D)[:-1]:
        du = dup_degree(R[i-1])
        dv = dup_degree(R[i  ])
        dw = dup_degree(R[i+1])

        if du % 2 and dv % 2:
            s = -s

        lc, i = dup_LC(R[i], K), i+1

        p *= b**dv * lc**(du-dw)
        q *= lc**(dv*(1+d))

    if s < 0:
        p = -p

    j = dup_degree(R[-2])

    res = dup_LC(R[-1], K)**j
    res = K.quo(res*p, q)

    return res, R

def dup_resultant(f, g, K):
    """Computes resultant of two polynomials in `K[x]`. """
    return dup_inner_resultant(f, g, K)[0]

def dmp_inner_subresultants(f, g, u, K):
    """Subresultant PRS algorithm in `K[X]`. """
    if not u:
        return dup_inner_subresultants(f, g, K)

    n = dmp_degree(f, u)
    m = dmp_degree(g, u)

    if n < m:
        f, g = g, f
        n, m = m, n

    R = [f, g]
    d = n - m
    v = u - 1

    b = dmp_pow(dmp_ground(-K.one, v), d+1, v, K)
    c = dmp_ground(-K.one, v)

    B, D = [b], [d]

    if dmp_zero_p(f, u) or dmp_zero_p(g, u):
        return R, B, D

    h = dmp_prem(f, g, u, K)
    h = dmp_mul_term(h, b, 0, u, K)

    while not dmp_zero_p(h, u):
        k = dmp_degree(h, u)
        R.append(h)

        lc = dmp_LC(g, K)

        p = dmp_pow(dmp_neg(lc, v, K), d, v, K)

        if not d:
            q = c
        else:
            q = dmp_pow(c, d-1, v, K)

        c = dmp_exquo(p, q, v, K)
        b = dmp_mul(dmp_neg(lc, v, K),
                    dmp_pow(c, m-k, v, K), v, K)

        f, g, m, d = g, h, k, m-k

        B.append(b)
        D.append(d)

        h = dmp_prem(f, g, u, K)
        h = [ dmp_exquo(ch, b, v, K) for ch in h ]

    return R, B, D

def dmp_subresultants(f, g, u, K):
    """Computes subresultant PRS of two polynomials in `K[X]`. """
    return dmp_inner_subresultants(f, g, u, K)[0]

def dmp_inner_resultant(f, g, u, K):
    """Resultant algorithm in `K[X]` using subresultant PRS. """
    if not u:
        return dup_inner_resultant(f, g, K)

    if dmp_zero_p(f, u) or dmp_zero_p(g, u):
        return (dmp_zero(u-1), [])

    R, B, D = dmp_inner_subresultants(f, g, u, K)

    if dmp_degree(R[-1], u) > 0:
        return (dmp_zero(u-1), R)
    if dmp_one_p(R[-2], u, K):
        return (dmp_LC(R[-1], K), R)

    s, i, v = 1, 1, u-1

    p = dmp_one(v, K)
    q = dmp_one(v, K)

    for b, d in zip(B, D)[:-1]:
        du = dmp_degree(R[i-1], u)
        dv = dmp_degree(R[i  ], u)
        dw = dmp_degree(R[i+1], u)

        if du % 2 and dv % 2:
            s = -s

        lc, i = dmp_LC(R[i], K), i+1

        p = dmp_mul(dmp_mul(p, dmp_pow(b, dv, v, K), v, K),
                               dmp_pow(lc, du-dw, v, K), v, K)
        q = dmp_mul(q, dmp_pow(lc, dv*(1+d), v, K), v, K)

        _, p, q = dmp_inner_gcd(p, q, v, K)

    if s < 0:
        p = dmp_neg(p, v, K)

    j = dmp_degree(R[-2], u)

    res = dmp_pow(dmp_LC(R[-1], K), j, v, K)
    res = dmp_quo(dmp_mul(res, p, v, K), q, v, K)

    return res, R

def dmp_resultant(f, g, u, K):
    """Computes resultant of two polynomials in `K[X]`. """
    return dmp_inner_resultant(f, g, u, K)[0]

def dup_discriminant(f, K):
    """Computes discriminant of a polynomial in `K[x]`. """
    d = dup_degree(f)

    if d <= 0:
        return K.zero
    else:
        s = (-1)**((d*(d-1)) // 2)
        c = dup_LC(f, K)

        r = dup_resultant(f, dup_diff(f, 1, K), K)

        return K.quo(r, c*s)

def dmp_discriminant(f, u, K):
    """Computes discriminant of a polynomial in `K[X]`. """
    if not u:
        return dup_discriminant(f, K)

    d, v = dmp_degree(f, u), u-1

    if d <= 0:
        return dmp_zero(v)
    else:
        s = (-1)**((d*(d-1)) // 2)
        c = dmp_LC(f, K)

        r = dmp_resultant(f, dmp_diff(f, 1, u, K), u, K)
        c = dmp_mul_ground(c, s, v, K)

        return dmp_quo(r, c, v, K)

def _dup_rr_trivial_gcd(f, g, K):
    """Handle trivial cases in GCD algorithm over a ring. """
    if not (f or g):
        return [], [], []
    elif not f:
        if K.is_nonnegative(dup_LC(g, K)):
            return g, [], [K.one]
        else:
            return dup_neg(g, K), [], [-K.one]
    elif not g:
        if K.is_nonnegative(dup_LC(f, K)):
            return f, [K.one], []
        else:
            return dup_neg(f, K), [-K.one], []

    return None

def _dup_ff_trivial_gcd(f, g, K):
    """Handle trivial cases in GCD algorithm over a field. """
    if not (f or g):
        return [], [], []
    elif not f:
        return dup_monic(g, K), [], [dup_LC(g, K)]
    elif not g:
        return dup_monic(f, K), [dup_LC(f, K)], []
    else:
        return None

USE_DMP_SIMPLIFY_GCD = 1

def _dmp_rr_trivial_gcd(f, g, u, K):
    """Handle trivial cases in GCD algorithm over a ring. """
    zero_f = dmp_zero_p(f, u)
    zero_g = dmp_zero_p(g, u)

    if zero_f and zero_g:
        return tuple(dmp_zeros(3, u, K))
    elif zero_f:
        if K.is_nonnegative(dmp_ground_LC(g, u, K)):
            return g, dmp_zero(u), dmp_one(u, K)
        else:
            return dmp_neg(g, u, K), dmp_zero(u), dmp_ground(-K.one, u)
    elif zero_g:
        if K.is_nonnegative(dmp_ground_LC(f, u, K)):
            return f, dmp_one(u, K), dmp_zero(u)
        else:
            return dmp_neg(f, u, K), dmp_ground(-K.one, u), dmp_zero(u)
    elif USE_DMP_SIMPLIFY_GCD:
        return _dmp_simplify_gcd(f, g, u, K)
    else:
        return None

def _dmp_ff_trivial_gcd(f, g, u, K):
    """Handle trivial cases in GCD algorithm over a field. """
    zero_f = dmp_zero_p(f, u)
    zero_g = dmp_zero_p(g, u)

    if zero_f and zero_g:
        return tuple(dmp_zeros(3, u, K))
    elif zero_f:
        return (dmp_ground_monic(g, u, K),
                dmp_zero(u),
                dmp_ground(dmp_ground_LC(g, u, K), u))
    elif zero_g:
        return (dmp_ground_monic(f, u, K),
                dmp_ground(dmp_ground_LC(f, u, K), u),
                dmp_zero(u))
    elif USE_DMP_SIMPLIFY_GCD:
        return _dmp_simplify_gcd(f, g, u, K)
    else:
        return None

def _dmp_simplify_gcd(f, g, u, K):
    """Try to eliminate `x_0` from GCD computation in `K[X]`. """
    df = dmp_degree(f, u)
    dg = dmp_degree(g, u)

    if df > 0 and dg > 0:
        return None

    if not (df or dg):
        F = dmp_LC(f, K)
        G = dmp_LC(g, K)
    else:
        if not df:
            F = dmp_LC(f, K)
            G = dmp_content(g, u, K)
        else:
            F = dmp_content(f, u, K)
            G = dmp_LC(g, K)

    v = u - 1
    h = dmp_gcd(F, G, v, K)

    cff = [ dmp_exquo(cf, h, v, K) for cf in f ]
    cfg = [ dmp_exquo(cg, h, v, K) for cg in g ]

    return [h], cff, cfg

def dup_rr_prs_gcd(f, g, K):
    """Computes polynomial GCD using subresultants over a ring. """
    result = _dup_rr_trivial_gcd(f, g, K)

    if result is not None:
        return result

    fc, F = dup_primitive(f, K)
    gc, G = dup_primitive(g, K)

    c = K.gcd(fc, gc)

    h = dup_subresultants(F, G, K)[-1]
    _, h = dup_primitive(h, K)

    if K.is_negative(dup_LC(h, K)):
        c = -c

    h = dup_mul_ground(h, c, K)

    cff = dup_exquo(f, h, K)
    cfg = dup_exquo(g, h, K)

    return h, cff, cfg

def dup_ff_prs_gcd(f, g, K):
    """Computes polynomial GCD using subresultants over a field. """
    result = _dup_ff_trivial_gcd(f, g, K)

    if result is not None:
        return result

    h = dup_subresultants(f, g, K)[-1]
    h = dup_monic(h, K)

    cff = dup_exquo(f, h, K)
    cfg = dup_exquo(g, h, K)

    return h, cff, cfg

def dmp_rr_prs_gcd(f, g, u, K):
    """Computes polynomial GCD using subresultants over a ring. """
    if not u:
        return dup_rr_prs_gcd(f, g, K)

    result = _dmp_rr_trivial_gcd(f, g, u, K)

    if result is not None:
        return result

    fc, F = dmp_primitive(f, u, K)
    gc, G = dmp_primitive(g, u, K)

    h = dmp_subresultants(F, G, u, K)[-1]
    c, _, _ = dmp_rr_prs_gcd(fc, gc, u-1, K)

    if K.is_negative(dmp_ground_LC(h, u, K)):
        h = dmp_neg(h, u, K)

    _, h = dmp_primitive(h, u, K)
    h = dmp_mul_term(h, c, 0, u, K)

    cff = dmp_exquo(f, h, u, K)
    cfg = dmp_exquo(g, h, u, K)

    return h, cff, cfg

def dmp_ff_prs_gcd(f, g, u, K):
    """Computes polynomial GCD using subresultants over a field. """
    if not u:
        return dup_ff_prs_gcd(f, g, K)

    result = _dmp_ff_trivial_gcd(f, g, u, K)

    if result is not None:
        return result

    fc, f = dmp_primitive(f, u, K)
    gc, g = dmp_primitive(g, u, K)

    h = dmp_subresultants(f, g, u, K)[-1]
    c, _, _ = dmp_ff_prs_gcd(fc, gc, u-1, K)

    _, h = dmp_primitive(h, u, K)
    h = dmp_mul_term(h, c, 0, u, K)
    h = dmp_ground_monic(h, u, K)

    cff = dmp_exquo(f, h, u, K)
    cfg = dmp_exquo(g, h, u, K)

    return h, cff, cfg

HEU_GCD_MAX = 6

def dup_zz_heu_gcd(f, g, K):
    """Heuristic polynomial GCD in `Z[x]`.

       Given univariate polynomials `f` and `g` in `Z[x]`, returns their GCD
       and cofactors, i.e. polynomials `h`, `cff` and `cfg` such that::

              h = gcd(f, g), cff = quo(f, h) and cfg = quo(g, h)

       The algorithm is purely heuristic which means it may fail to compute
       the GCD. This will be signaled by raising an exception. In this case
       you will need to switch to another GCD method.

       The algorithm computes the polynomial GCD by evaluating polynomials
       f and g at certain points and computing (fast) integer GCD of those
       evaluations. The polynomial GCD is recovered from the integer image
       by interpolation.  The final step is to verify if the result is the
       correct GCD. This gives cofactors as a side effect.

       References
       ==========

       .. [Liao95] Hsin-Chao Liao,  R. Fateman, Evaluation of the heuristic
          polynomial GCD, International Symposium on Symbolic and Algebraic
          Computation (ISSAC), ACM Press, Montreal, Quebec, Canada, 1995,
          pp. 240--247

    """
    result = _dup_rr_trivial_gcd(f, g, K)

    if result is not None:
        return result

    def interpolate(h, x):
        f = []

        while h:
            g = h % x

            if g > x // 2:
                g -= x

            f.insert(0, g)
            h = (h-g) // x

        return f

    def finalize(h, cff, cfg, gcd):
        h = dup_mul_ground(h, gcd, K)
        return h, cff, cfg

    df = dup_degree(f)
    dg = dup_degree(g)

    gcd, f, g = dup_extract(f, g, K)

    if df == 0 or dg == 0:
        return [gcd], f, g

    f_norm = dup_max_norm(f, K)
    g_norm = dup_max_norm(g, K)

    B = 2*min(f_norm, g_norm) + 29

    x = max(min(B, 99*K.sqrt(B)),
            2*min(f_norm // abs(dup_LC(f, K)),
                  g_norm // abs(dup_LC(g, K))) + 2)

    for i in xrange(0, HEU_GCD_MAX):
        ff = dup_eval(f, x, K)
        gg = dup_eval(g, x, K)

        if ff and gg:
            h = K.gcd(ff, gg)

            cff = ff // h
            cfg = gg // h

            h = interpolate(h, x)
            h = dup_primitive(h, K)[1]

            cff_, r = dup_div(f, h, K)

            if not r:
                cfg_, r = dup_div(g, h, K)

                if not r:
                    return finalize(h, cff_, cfg_, gcd)

            cff = interpolate(cff, x)

            h, r = dup_div(f, cff, K)

            if not r:
                cfg_, r = dup_div(g, h, K)

                if not r:
                    return finalize(h, cff, cfg_, gcd)

            cfg = interpolate(cfg, x)

            h, r = dup_div(g, cfg, K)

            if not r:
                cff_, r = dup_div(f, h, K)

                if not r:
                    return finalize(h, cff_, cfg, gcd)

        x = 73794*x * K.sqrt(K.sqrt(x)) // 27011

    raise HeuristicGCDFailed('no luck')

def dmp_zz_heu_gcd(f, g, u, K):
    """Heuristic polynomial GCD in `Z[X]`.

       Given univariate polynomials `f` and `g` in `Z[X]`, returns their GCD
       and cofactors, i.e. polynomials `h`, `cff` and `cfg` such that::

              h = gcd(f, g), cff = quo(f, h) and cfg = quo(g, h)

       The algorithm is purely heuristic which means it may fail to compute
       the GCD. This will be signaled by raising an exception. In this case
       you will need to switch to another GCD method.

       The algorithm computes the polynomial GCD by evaluating polynomials
       f and g at certain points and computing (fast) integer GCD of those
       evaluations. The polynomial GCD is recovered from the integer image
       by interpolation. The evaluation proces reduces f and g variable by
       variable into a large integer.  The final step  is to verify if the
       interpolated polynomial is the correct GCD. This gives cofactors of
       the input polynomials as a side effect.

       References
       ==========

       .. [Liao95] Hsin-Chao Liao,  R. Fateman, Evaluation of the heuristic
          polynomial GCD, International Symposium on Symbolic and Algebraic
          Computation (ISSAC), ACM Press, Montreal, Quebec, Canada, 1995,
          pp. 240--247

    """
    if not u:
        return dup_zz_heu_gcd(f, g, K)

    result = _dmp_rr_trivial_gcd(f, g, u, K)

    if result is not None:
        return result

    def interpolate(h, x, v):
        f = []

        while not dmp_zero_p(h, v):
            g = dmp_ground_trunc(h, x, v, K)
            f.insert(0, g)

            h = dmp_sub(h, g, v, K)
            h = dmp_exquo_ground(h, x, v, K)

        if K.is_negative(dmp_ground_LC(f, v+1, K)):
            return dmp_neg(f, v+1, K)
        else:
            return f

    def finalize(h, cff, cfg, gcd):
        h = dmp_mul_ground(h, gcd, u, K)
        return h, cff, cfg

    df = dmp_degree(f, u)
    dg = dmp_degree(g, u)

    gcd, f, g = dmp_ground_extract(f, g, u, K)

    f_norm = dmp_max_norm(f, u, K)
    g_norm = dmp_max_norm(g, u, K)

    B = 2*min(f_norm, g_norm) + 29

    x = max(min(B, 99*K.sqrt(B)),
            2*min(f_norm // abs(dmp_ground_LC(f, u, K)),
                  g_norm // abs(dmp_ground_LC(g, u, K))) + 2)

    for i in xrange(0, HEU_GCD_MAX):
        ff = dmp_eval(f, x, u, K)
        gg = dmp_eval(g, x, u, K)

        v = u - 1

        if not (dmp_zero_p(ff, v) or dmp_zero_p(gg, v)):
            h, cff, cfg = dmp_zz_heu_gcd(ff, gg, v, K)

            h = interpolate(h, x, v)
            h = dmp_ground_primitive(h, u, K)[1]

            cff_, r = dmp_div(f, h, u, K)

            if dmp_zero_p(r, u):
                cfg_, r = dmp_div(g, h, u, K)

                if dmp_zero_p(r, u):
                    return finalize(h, cff_, cfg_, gcd)

            cff = interpolate(cff, x, v)

            h, r = dmp_div(f, cff, u, K)

            if dmp_zero_p(r, u):
                cfg_, r = dmp_div(g, h, u, K)

                if dmp_zero_p(r, u):
                    return finalize(h, cff, cfg_, gcd)

            cfg = interpolate(cfg, x, v)

            h, r = dmp_div(g, cfg, u, K)

            if dmp_zero_p(r, u):
                cff_, r = dmp_div(f, h, u, K)

                if dmp_zero_p(r, u):
                    return finalize(h, cff_, cfg, gcd)

        x = 73794*x * K.sqrt(K.sqrt(x)) // 27011

    raise HeuristicGCDFailed('no luck')

def dup_qq_heu_gcd(f, g, K0):
    """Heuristic polynomial GCD in `Q[x]`. """
    result = _dup_ff_trivial_gcd(f, g, K0)

    if result is not None:
        return result

    K1 = K0.get_ring()

    cf, f = dup_ground_to_ring(f, K0, K1)
    cg, g = dup_ground_to_ring(g, K0, K1)

    f = dup_convert(f, K0, K1)
    g = dup_convert(g, K0, K1)

    h, cff, cfg = dup_zz_heu_gcd(f, g, K1)

    h = dup_convert(h, K1, K0)

    c = dup_LC(h, K0)
    h = dup_monic(h, K0)

    cff = dup_convert(cff, K1, K0)
    cfg = dup_convert(cfg, K1, K0)

    cff = dup_mul_ground(cff, K0.quo(c, cf), K0)
    cfg = dup_mul_ground(cfg, K0.quo(c, cg), K0)

    return h, cff, cfg

def dmp_qq_heu_gcd(f, g, u, K0):
    """Heuristic polynomial GCD in `Q[X]`. """
    result = _dmp_ff_trivial_gcd(f, g, u, K0)

    if result is not None:
        return result

    K1 = K0.get_ring()

    cf, f = dmp_ground_to_ring(f, u, K0, K1)
    cg, g = dmp_ground_to_ring(g, u, K0, K1)

    f = dmp_convert(f, u, K0, K1)
    g = dmp_convert(g, u, K0, K1)

    h, cff, cfg = dmp_zz_heu_gcd(f, g, u, K1)

    h = dmp_convert(h, u, K1, K0)

    c = dmp_ground_LC(h, u, K0)
    h = dmp_ground_monic(h, u, K0)

    cff = dmp_convert(cff, u, K1, K0)
    cfg = dmp_convert(cfg, u, K1, K0)

    cff = dmp_mul_ground(cff, K0.quo(c, cf), u, K0)
    cfg = dmp_mul_ground(cfg, K0.quo(c, cg), u, K0)

    return h, cff, cfg

USE_DUP_HEU_GCD = 1
USE_DMP_HEU_GCD = 1

def dup_inner_gcd(f, g, K):
    """Computes polynomial GCD and cofactors of `f` and `g` in `K[x]`. """
    if K.has_Field:
        if USE_DUP_HEU_GCD:
            if K.is_QQ:
                try:
                    return dup_qq_heu_gcd(f, g, K)
                except HeuristicGCDFailed:
                    pass

        return dup_ff_prs_gcd(f, g, K)
    else:
        if USE_DUP_HEU_GCD:
            if K.is_ZZ:
                try:
                    return dup_zz_heu_gcd(f, g, K)
                except HeuristicGCDFailed:
                    pass

        return dup_rr_prs_gcd(f, g, K)

def _dmp_inner_gcd(f, g, u, K):
    """Helper function for `dmp_inner_gcd()`. """
    if K.has_Field:
        if USE_DMP_HEU_GCD:
            if K.is_QQ:
                try:
                    return dmp_qq_heu_gcd(f, g, u, K)
                except HeuristicGCDFailed:
                    pass

        return dmp_ff_prs_gcd(f, g, u, K)
    else:
        if USE_DMP_HEU_GCD:
            if K.is_ZZ:
                try:
                     return dmp_zz_heu_gcd(f, g, u, K)
                except HeuristicGCDFailed:
                    pass

        return dmp_rr_prs_gcd(f, g, u, K)

def dmp_inner_gcd(f, g, u, K):
    """Computes polynomial GCD and cofactors of `f` and `g` in `K[X]`. """
    if not u:
        return dup_inner_gcd(f, g, K)

    J, (f, g) = dmp_multi_deflate((f, g), u, K)
    h, cff, cfg = _dmp_inner_gcd(f, g, u, K)

    return (dmp_inflate(h, J, u, K),
            dmp_inflate(cff, J, u, K),
            dmp_inflate(cfg, J, u, K))

def dup_gcd(f, g, K):
    """Computes polynomial GCD of `f` and `g` in `K[x]`. """
    return dup_inner_gcd(f, g, K)[0]

def dmp_gcd(f, g, u, K):
    """Computes polynomial GCD of `f` and `g` in `K[X]`. """
    return dmp_inner_gcd(f, g, u, K)[0]

def dup_rr_lcm(f, g, K):
    """Computes polynomial LCM over a ring in `K[x]`. """
    fc, f = dup_primitive(f, K)
    gc, g = dup_primitive(g, K)

    c = K.lcm(fc, gc)

    h = dup_exquo(dup_mul(f, g, K),
                  dup_gcd(f, g, K), K)

    return dup_mul_ground(h, c, K)

def dup_ff_lcm(f, g, K):
    """Computes polynomial LCM over a field in `K[x]`. """
    h = dup_exquo(dup_mul(f, g, K),
                  dup_gcd(f, g, K), K)

    return dup_ground_monic(h, K)

def dup_lcm(f, g, K):
    """Computes polynomial LCM of `f` and `g` in `K[x]`. """
    if K.has_Field:
        return dup_ff_lcm(f, g, K)
    else:
        return dup_rr_lcm(f, g, K)

def dmp_rr_lcm(f, g, u, K):
    """Computes polynomial LCM over a ring in `K[X]`. """
    fc, f = dmp_ground_primitive(f, u, K)
    gc, g = dmp_ground_primitive(g, u, K)

    c = K.lcm(fc, gc)

    h = dmp_exquo(dmp_mul(f, g, u, K),
                  dmp_gcd(f, g, u, K), u, K)

    return dmp_mul_ground(h, c, u, K)

def dmp_ff_lcm(f, g, u, K):
    """Computes polynomial LCM over a field in `K[X]`. """
    h = dmp_exquo(dmp_mul(f, g, u, K),
                  dmp_gcd(f, g, u, K), u, K)

    return dmp_ground_monic(h, u, K)

def dmp_lcm(f, g, u, K):
    """Computes polynomial LCM of `f` and `g` in `K[X]`. """
    if not u:
        return dup_lcm(f, g, K)

    if K.has_Field:
        return dmp_ff_lcm(f, g, u, K)
    else:
        return dmp_rr_lcm(f, g, u, K)

def dup_trunc(f, p, K):
    """Reduce `K[x]` polynomial modulo a constant `p` in `K`. """
    if K.is_ZZ:
        g = []

        for c in f:
            c = c % p

            if c > p // 2:
                g.append(c - p)
            else:
                g.append(c)
    else:
        g = [ c % p for c in f ]

    return dup_strip(g)

def dmp_trunc(f, p, u, K):
    """Reduce `K[X]` polynomial modulo a polynomial `p` in `K[Y]`. """
    return dmp_strip([ dmp_rem(c, p, u-1, K) for c in f ], u)

def dmp_ground_trunc(f, p, u, K):
    """Reduce `K[X]` polynomial modulo a constant `p` in `K`. """
    if not u:
        return dup_trunc(f, p, K)
    else:
        return dmp_strip([ dmp_ground_trunc(c, p, u-1, K) for c in f ], u)

def dup_monic(f, K):
    """Divides all coefficients by `LC(f)` in `K[x]`. """
    if not f:
        return f

    lc = dup_LC(f, K)

    if K.is_one(lc):
        return f
    else:
        return dup_quo_ground(f, lc, K)

def dmp_ground_monic(f, u, K):
    """Divides all coefficients by `LC(f)` in `K[X]`. """
    if not u:
        return dup_monic(f, K)

    if dmp_zero_p(f, u):
        return f

    lc = dmp_ground_LC(f, u, K)

    if K.is_one(lc):
        return f
    else:
        return dmp_quo_ground(f, lc, u, K)

def dup_rr_content(f, K):
    """Returns GCD of coefficients over a ring. """
    cont = K.zero

    for c in f:
        cont = K.gcd(cont, c)

        if K.is_one(cont):
            break

    return cont

def dup_ff_content(f, K):
    """Returns GCD of coefficients over a field. """
    if not f:
        return K.zero
    else:
        return K.one

def dup_content(f, K):
    """Returns GCD of coefficients in `K[x]`. """
    if K.has_Field:
        return dup_ff_content(f, K)
    else:
        return dup_rr_content(f, K)

def dmp_content(f, u, K):
    """Returns GCD of multivariate coefficients. """
    cont, v = dmp_LC(f, K), u-1

    if dmp_zero_p(f, u):
        return cont

    for c in f[1:]:
        cont = dmp_gcd(cont, c, v, K)

        if dmp_one_p(cont, v, K):
            break

    if K.is_negative(dmp_ground_LC(cont, v, K)):
        return dmp_neg(cont, v, K)
    else:
        return cont

def dmp_rr_ground_content(f, u, K):
    """Returns GCD of coefficients over a ring. """
    if not u:
        return dup_rr_content(f, K)

    cont, v = K.zero, u-1

    for c in f:
        gc = dmp_rr_ground_content(c, v, K)
        cont = K.gcd(cont, gc)

        if K.is_one(cont):
            break

    return cont

def dmp_ff_ground_content(f, u, K):
    """Returns GCD of coefficients over a field. """
    if not f:
        return K.zero
    else:
        return K.one

def dmp_ground_content(f, u, K):
    """Returns GCD of coefficients in `K[X]`. """
    if not u:
        return dup_content(f, K)

    if K.has_Field:
        return dmp_ff_ground_content(f, u, K)
    else:
        return dmp_rr_ground_content(f, u, K)

def dup_rr_primitive(f, K):
    """Returns content and a primitive polynomial over a ring. """
    cont = dup_content(f, K)

    if not f or K.is_one(cont):
        return cont, f
    else:
        return cont, dup_exquo_ground(f, cont, K)

def dup_ff_primitive(f, K):
    """Returns content and a primitive polynomial over a field. """
    return K.one, f

def dup_primitive(f, K):
    """Returns content and a primitive polynomial in `K[x]`. """
    if K.has_Field:
        return dup_ff_primitive(f, K)
    else:
        return dup_rr_primitive(f, K)

def dmp_primitive(f, u, K):
    """Returns multivariate content and a primitive polynomial. """
    cont, v = dmp_content(f, u, K), u-1

    if dmp_zero_p(f, u) or dmp_one_p(cont, v, K):
        return cont, f
    else:
        return cont, [ dmp_exquo(c, cont, v, K) for c in f ]

def dmp_rr_ground_primitive(f, u, K):
    """Returns content and a primitive polynomial over a ring. """
    cont = dmp_ground_content(f, u, K)

    if K.is_one(cont):
        return cont, f
    else:
        return cont, dmp_exquo_ground(f, cont, u, K)

def dmp_ff_ground_primitive(f, u, K):
    """Returns content and a primitive polynomial over a ring. """
    if dmp_zero_p(f, u):
        return K.zero, f
    else:
        return K.one, f

def dmp_ground_primitive(f, u, K):
    """Returns content and a primitive polynomial in `K[x]`. """
    if not u:
        return dup_primitive(f, K)

    if dmp_zero_p(f, u):
        return K.zero, f

    if K.has_Field:
        return dmp_ff_ground_primitive(f, u, K)
    else:
        return dmp_rr_ground_primitive(f, u, K)

def dup_sqf_p(f, K):
    """Returns `True` if `f` is a square-free polynomial in `K[x]`. """
    if not f:
        return True
    else:
        return not dup_degree(dup_gcd(f, dup_diff(f, 1, K), K))

def dmp_sqf_p(f, u, K):
    """Returns `True` if `f` is a square-free polynomial in `K[X]`. """
    if dmp_zero_p(f, u):
        return True
    else:
        return not dmp_degree(dmp_gcd(f, dmp_diff(f, 1, u, K), u, K), u)

def dup_sqf_norm(f, K):
    """Square-free norm of `f`, useful over algebraic domains. """
    s, g = 0, dmp_raise(K.mod, 1, K.dom)

    while True:
        h, _ = dmp_inject(f, K, front=True)
        r = dmp_resultant(g, h, 1, K.dom)

        if dup_sqf_p(r, K.dom):
            return s, h, r
        else:
            f = dup_compose(f, [K.one,-K.alpha], K)

def dmp_sqf_norm(f, u, K):
    """Square-free norm of `f`, useful over algebraic domains. """
    raise NotImplementedError('algebraic numbers')

def dup_sqf_part(f, K):
    """Returns square-free part of a polynomial in `K[x]`. """
    if not f:
        return f

    if K.is_negative(dup_LC(f, K)):
        f = dup_neg(f, K)

    gcd = dup_gcd(f, dup_diff(f, 1, K), K)
    sqf = dup_exquo(f, gcd, K)

    if K.has_Field:
        return dup_monic(sqf, K)
    else:
        return dup_primitive(sqf, K)[1]

def dmp_sqf_part(f, u, K):
    """Returns square-free part of a polynomial in `K[X]`. """
    if dmp_zero_p(f, u):
        return f

    if K.is_negative(dmp_ground_LC(f, u, K)):
        f = dmp_neg(f, u, K)

    gcd = dmp_gcd(f, dmp_diff(f, 1, u, K), u, K)
    sqf = dmp_exquo(f, gcd, u, K)

    if K.has_Field:
        return dmp_ground_monic(sqf, u, K)
    else:
        return dmp_primitive(sqf, u, K)[1]

def dup_sqf_list(f, K, **args):
    """Returns square-free decomposition of a polynomial in `K[x]`. """
    if K.has_Field:
        coeff = dup_LC(f, K)
        f = dup_monic(f, K)
    else:
        coeff, f = dup_primitive(f, K)

        if K.is_negative(dup_LC(f, K)):
            f = dup_neg(f, K)
            coeff = -coeff

    if dup_degree(f) <= 0:
        if args.get('include', False):
            return f
        else:
            return coeff, []

    result, i = [], 1

    h = dup_diff(f, 1, K)
    g, p, q = dup_inner_gcd(f, h, K)

    all = args.get('all', False)

    while True:
        d = dup_diff(p, 1, K)
        h = dup_sub(q, d, K)

        if not h:
            result.append((p, i))
            break

        g, p, q = dup_inner_gcd(p, h, K)

        if all or dup_degree(g) > 0:
            result.append((g, i))

        i += 1

    if not args.get('include', False):
        return coeff, result
    else:
        (g, i), rest = result[0], result[1:]
        g = dup_mul_ground(g, coeff, K)

        return [(g, i)] + rest

def dmp_sqf_list(f, u, K, **args):
    """Returns square-free decomposition of a polynomial in `K[X]`. """
    if not u:
        return dup_sqf_list(f, K, **args)

    if K.has_Field:
        coeff = dmp_ground_LC(f, u, K)
        f = dmp_ground_monic(f, u, K)
    else:
        coeff, f = dmp_ground_primitive(f, u, K)

        if K.is_negative(dmp_ground_LC(f, u, K)):
            f = dmp_neg(f, u, K)
            coeff = -coeff

    if dmp_degree(f, u) <= 0:
        if args.get('include', False):
            return f
        else:
            return coeff, []

    result, i = [], 1

    h = dmp_diff(f, 1, u, K)
    g, p, q = dmp_inner_gcd(f, h, u, K)

    all = args.get('all', False)

    while True:
        d = dmp_diff(p, 1, u, K)
        h = dmp_sub(q, d, u, K)

        if dmp_zero_p(h, u):
            result.append((p, i))
            break

        g, p, q = dmp_inner_gcd(p, h, u, K)

        if all or dmp_degree(g, u) > 0:
            result.append((g, i))

        i += 1

    if not args.get('include', False):
        return coeff, result
    else:
        (g, i), rest = result[0], result[1:]
        g = dup_mul_ground(g, coeff, K)

        return [(g, i)] + rest

def dup_extract(f, g, K):
    """Extracts common content from a pair of polynomials in `K[x]`. """
    fc = dup_content(f, K)
    gc = dup_content(g, K)

    gcd = K.gcd(fc, gc)

    if not K.is_one(gcd):
        f = dup_exquo_ground(f, gcd, K)
        g = dup_exquo_ground(g, gcd, K)

    return gcd, f, g

def dmp_ground_extract(f, g, u, K):
    """Extracts common content from a pair of polynomials in `K[X]`. """
    fc = dmp_ground_content(f, u, K)
    gc = dmp_ground_content(g, u, K)

    gcd = K.gcd(fc, gc)

    if not K.is_one(gcd):
        f = dmp_exquo_ground(f, gcd, u, K)
        g = dmp_exquo_ground(g, gcd, u, K)

    return gcd, f, g

def dup_compose(f, g, K):
    """Evaluate functional composition `f(g)` in `K[x]`. """
    if len(g) <= 1:
        return dup_strip([dup_eval(f, dup_LC(g, K), K)])

    if not f:
        return []

    h = [f[0]]

    for c in f[1:]:
        h = dup_mul(h, g, K)
        h = dup_add_term(h, c, 0, K)

    return h

def dup_decompose(f, K):
    """Computes functional decomposition of `f` in `K[x]`.

       Given an univariate polynomial `f` with coefficients in a field of
       characteristic zero, returns tuple `(f_1, f_2, ..., f_n)`, where::

                  f = f_1 o f_2 o ... f_n = f_1(f_2(... f_n))

       and `f_2, ..., f_n` are monic and homogeneous polynomials of at
       least second degree.

       Unlike factorization, complete functional decompositions of
       polynomials are not unique, consider examples:

       1. `f o g = f(x + b) o (g - b)`
       2. `x**n o x**m = x**m o x**n`
       3. `T_n o T_m = T_m o T_n`

       where `T_n` and `T_m` are Chebyshev polynomials.

       References
       ==========

       .. [Kozen89] D. Kozen, S. Landau, Polynomial decomposition algorithms,
          Journal of Symbolic Computation 7 (1989), pp. 445-456

    """
    def right_decompose(f, s):
        n = dup_degree(f)
        lc = dup_LC(f, K)

        f = dup_to_raw_dict(f)
        g = { s : K.one }

        r = n // s

        for i in xrange(1, s):
            coeff = K.zero

            for j in xrange(0, i):
                if not n+j-i in f:
                    continue

                if not s-j in g:
                    continue

                fc, gc = f[n+j-i], g[s-j]
                coeff += (i - r*j)*fc*gc

            g[s-i] = K.quo(coeff, i*r*lc)

        return dup_from_raw_dict(g, K)

    def left_decompose(f, h):
        g, i = {}, 0

        while f:
            q, r = dup_div(f, h, K)

            if dup_degree(r) > 0:
                return None
            else:
                g[i] = dup_LC(r, K)
                f, i = q, i + 1

        return dup_from_raw_dict(g, K)

    def decompose(f):
        df = dup_degree(f)

        for s in xrange(2, df):
            if df % s != 0:
                continue

            h = right_decompose(f, s)

            if h is not None:
                g = left_decompose(f, h)

                if g is not None:
                    return g, h

        return None

    F = []

    while True:
        result = decompose(f)

        if result is not None:
            f, h = result
            F = [h] + F
        else:
            break

    return [f] + F

def dup_sturm(f, K):
    """Computes the Sturm sequence of `f` in `F[x]`.

       Given an univariate, square-free polynomial `f(x)` returns the
       associated Sturm sequence `f_0(x), ..., f_n(x)` defined by::

           f_0(x), f_1(x) = f(x), f'(x)
           f_n = -rem(f_{n-2}(x), f_{n-1}(x))

       References
       ==========

       .. [Davenport88] J.H. Davenport, Y. Siret, E. Tournier,
           Computer Algebra Systems and Algorithms for Algebraic
           Computation, Academic Press, London, 1988, pp. 124-128

    """
    if not K.has_Field:
        raise DomainError('computation can be done only in a field')

    f = dup_sqf_part(f, K)

    sturm = [f, dup_diff(f, 1, K)]

    while sturm[-1]:
        s = dup_rem(sturm[-2], sturm[-1], K)
        sturm.append(dup_neg(s, K))

    return sturm[:-1]

