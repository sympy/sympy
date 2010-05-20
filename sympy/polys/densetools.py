"""Advanced tools for dense recursive polynomials in `K[x]` or `K[X]`. """

from sympy.polys.densebasic import (
    dup_strip, dmp_strip,
    dup_reverse,
    dup_convert, dmp_convert,
    dup_degree, dmp_degree, dmp_degree_in,
    dup_to_dict, dmp_to_dict,
    dup_from_dict, dmp_from_dict,
    dup_LC, dmp_LC, dmp_ground_LC,
    dup_TC, dmp_TC, dmp_ground_TC,
    dmp_zero, dmp_one, dmp_ground,
    dmp_zero_p, dmp_one_p,
    dmp_multi_deflate, dmp_inflate,
    dup_to_raw_dict, dup_from_raw_dict,
    dmp_raise, dmp_apply_pairs,
    dmp_inject, dmp_zeros,
    dup_terms_gcd
)

from sympy.polys.densearith import (
    dup_add_term, dmp_add_term,
    dup_mul_term, dmp_mul_term,
    dup_lshift, dup_rshift,
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
    dup_expand, dmp_expand,
    dup_add_mul, dup_sub_mul,
    dup_mul_ground, dmp_mul_ground,
    dup_quo_ground, dmp_quo_ground,
    dup_exquo_ground, dmp_exquo_ground,
    dup_max_norm, dmp_max_norm
)

from sympy.polys.galoistools import (
    gf_int, gf_crt
)

from sympy.polys.polyerrors import (
    HeuristicGCDFailed,
    HomomorphismFailed,
    RefinementFailed,
    NotInvertible,
    DomainError
)

from sympy.polys.polyconfig import query

from sympy.ntheory import nextprime

from sympy.utilities import (
    cythonized, variations
)

from random import random as randfloat
from operator import itemgetter

def dup_clear_denoms(f, K0, K1=None, convert=False):
    """Clear denominators, i.e. transform `K_0` to `K_1`. """
    if K1 is None:
        K1 = K0.get_ring()

    common = K1.one

    for c in f:
        common = K1.lcm(common, K0.denom(c))

    if not K1.is_one(common):
        f = dup_mul_ground(f, common, K0)

    if not convert:
        return common, f
    else:
        return common, dup_convert(f, K0, K1)

@cythonized("v,w")
def _rec_clear_denoms(g, v, K0, K1):
    """XXX"""
    common = K1.one

    if not v:
        for c in g:
            common = K1.lcm(common, K0.denom(c))
    else:
        w = v-1

        for c in g:
            common = K1.lcm(common, _rec_clear_denoms(c, w, K0, K1))

    return common

@cythonized("u")
def dmp_clear_denoms(f, u, K0, K1=None, convert=False):
    """Clear denominators, i.e. transform `K_0` to `K_1`. """
    if not u:
        return dup_clear_denoms(f, K0, K1)

    if K1 is None:
        K1 = K0.get_ring()

    common = _rec_clear_denoms(f, u, K0, K1)

    if not K1.is_one(common):
        f = dmp_mul_ground(f, common, u, K0)

    if not convert:
        return common, f
    else:
        return common, dmp_convert(f, u, K0, K1)

@cythonized("m,n,i,j")
def dup_integrate(f, m, K):
    """Computes indefinite integral of `f` in `K[x]`. """
    if m <= 0 or not f:
        return f

    g = [K.zero]*m

    for i, c in enumerate(reversed(f)):
        n = i+1

        for j in xrange(1, m):
            n *= i+j+1

        g.insert(0, K.quo(c, K(n)))

    return g

@cythonized("m,u,v,n,i,j")
def dmp_integrate(f, m, u, K):
    """Computes indefinite integral of `f` in `x_0` in `K[X]`. """
    if not u:
        return dup_integrate(f, m, K)

    if m <= 0 or dmp_zero_p(f, u):
        return f

    g, v = dmp_zeros(m, u-1, K), u-1

    for i, c in enumerate(reversed(f)):
        n = i+1

        for j in xrange(1, m):
            n *= i+j+1

        g.insert(0, dmp_quo_ground(c, K(n), v, K))

    return g

@cythonized("m,v,w,i,j")
def _rec_integrate_in(g, m, v, i, j, K):
    """XXX"""
    if i == j:
        return dmp_integrate(g, m, v, K)

    w, i = v-1, i+1

    return dmp_strip([ _rec_integrate_in(c, m, w, i, j, K) for c in g ], v)

@cythonized("m,j,u")
def dmp_integrate_in(f, m, j, u, K):
    """Computes indefinite integral of `f` in `x_j` in `K[X]`. """
    if j < 0 or j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))

    return _rec_integrate_in(f, m, u, 0, j, K)

@cythonized("m,n,i")
def dup_diff(f, m, K):
    """m-th order derivative of a polynomial in `K[x]`. """
    if m <= 0:
        return f

    n = dup_degree(f)

    if n < m:
        return []

    deriv, c = [], K.one

    for i in xrange(0, m):
        c, n = c*K(n), n-1

    for coeff in f[:-m]:
        deriv.append(coeff*c)
        c, n = K(n)*K.exquo(c, K(n+m)), n-1

    return deriv

@cythonized("u,v,m,n,i")
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
        c, n = c*K(n), n-1

    for coeff in f[:-m]:
        h = dmp_mul_ground(coeff, c, v, K)
        c, n = K(n)*K.exquo(c, K(n+m)), n-1
        deriv.append(h)

    return deriv

@cythonized("m,v,w,i,j")
def _rec_diff_in(g, m, v, i, j, K):
    """XXX"""
    if i == j:
        return dmp_diff(g, m, v, K)

    w, i = v-1, i+1

    return dmp_strip([ _rec_diff_in(c, m, w, i, j, K) for c in g ], v)

@cythonized("m,j,u")
def dmp_diff_in(f, m, j, u, K):
    """m-th order derivative in `x_j` of a polynomial in `K[X]`. """
    if j < 0 or j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))

    return _rec_diff_in(f, m, u, 0, j, K)

def dup_eval(f, a, K):
    """Evaluate a polynomial at `x = a` in `K[x]` using Horner scheme. """
    if not a:
        return dup_TC(f, K)

    result = K.zero

    for c in f:
        result *= a
        result += c

    return result

@cythonized("u,v")
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

@cythonized("v,i,j")
def _rec_eval_in(g, a, v, i, j, K):
    """XXX"""
    if i == j:
        return dmp_eval(g, a, v, K)

    v, i = v-1, i+1

    return dmp_strip([ _rec_eval_in(c, a, v, i, j, K) for c in g ], v)

@cythonized("u")
def dmp_eval_in(f, a, j, u, K):
    """Evaluate a polynomial at `x_j = a` in `K[X]` using Horner scheme. """
    if j < 0 or j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))

    return _rec_eval_in(f, a, u, 0, j, K)

@cythonized("i,u")
def _rec_eval_tail(g, i, A, u, K):
    """XXX"""
    if i == u:
        return dup_eval(g, A[-1], K)
    else:
        h = [ _rec_eval_tail(c, i+1, A, u, K) for c in g ]

        if i < u - len(A) + 1:
            return h
        else:
            return dup_eval(h, A[-u+i-1], K)

@cythonized("u")
def dmp_eval_tail(f, A, u, K):
    """Evaluate a polynomial at `x_j = a_j, ...` in `K[X]`. """
    if not A:
        return f

    if dmp_zero_p(f, u):
        return dmp_zero(u - len(A))

    e = _rec_eval_tail(f, 0, A, u, K)

    if u == len(A)-1:
        return e
    else:
        return dmp_strip(e, u - len(A))

@cythonized("m,v,i,j")
def _rec_diff_eval(g, m, a, v, i, j, K):
    """XXX"""
    if i == j:
        return dmp_eval(dmp_diff(g, m, v, K), a, v, K)

    v, i = v-1, i+1

    return dmp_strip([ _rec_diff_eval(c, m, a, v, i, j, K) for c in g ], v)

@cythonized("m,j,u")
def dmp_diff_eval_in(f, m, a, j, u, K):
    """Differentiate and evaluate a polynomial in `x_j` at `a` in `K[X]`. """
    if j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))
    if not j:
        return dmp_eval(dmp_diff(f, m, u, K), a, u, K)

    return _rec_diff_eval(f, m, a, u, 0, j, K)

def dup_half_gcdex(f, g, K):
    """Half extended Euclidean algorithm in `F[x]`. """
    if not (K.has_Field or not K.is_Exact):
        raise DomainError("can't compute half extended GCD over %s" % K)

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

@cythonized("n,m,d,k")
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

@cythonized("s,i,du,dv,dw")
def dup_prs_resultant(f, g, K):
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

    i = dup_degree(R[-2])

    res = dup_LC(R[-1], K)**i
    res = K.quo(res*p, q)

    return res, R

def dup_resultant(f, g, K):
    """Computes resultant of two polynomials in `K[x]`. """
    return dup_prs_resultant(f, g, K)[0]

@cythonized("u,v,n,m,d,k")
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

@cythonized("u")
def dmp_subresultants(f, g, u, K):
    """Computes subresultant PRS of two polynomials in `K[X]`. """
    return dmp_inner_subresultants(f, g, u, K)[0]

@cythonized("u,v,s,i,d,du,dv,dw")
def dmp_prs_resultant(f, g, u, K):
    """Resultant algorithm in `K[X]` using subresultant PRS. """
    if not u:
        return dup_prs_resultant(f, g, K)

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

    i = dmp_degree(R[-2], u)

    res = dmp_pow(dmp_LC(R[-1], K), i, v, K)
    res = dmp_quo(dmp_mul(res, p, v, K), q, v, K)

    return res, R

@cythonized("u,v,n,m,N,M,B")
def dmp_zz_modular_resultant(f, g, p, u, K):
    """Compute resultant of `f` and `g` modulo a prime `p`. """
    if not u:
        return gf_int(dup_prs_resultant(f, g, K)[0] % p, p)

    v = u - 1

    n = dmp_degree(f, u)
    m = dmp_degree(g, u)

    N = dmp_degree_in(f, 1, u)
    M = dmp_degree_in(g, 1, u)

    B = n*M + m*N

    D, a = [K.one], -K.one
    r = dmp_zero(v)

    while dup_degree(D) <= B:
        while True:
            a += K.one

            if a == p:
                raise HomomorphismFailed('no luck')

            F = dmp_eval_in(f, gf_int(a, p), 1, u, K)

            if dmp_degree(F, v) == n:
                G = dmp_eval_in(g, gf_int(a, p), 1, u, K)

                if dmp_degree(G, v) == m:
                    break

        R = dmp_zz_modular_resultant(F, G, p, v, K)
        e = dmp_eval(r, a, v, K)

        if not v:
            R = dup_strip([R])
            e = dup_strip([e])
        else:
            R = [R]
            e = [e]

        d = K.invert(dup_eval(D, a, K), p)
        d = dup_mul_ground(D, d, K)
        d = dmp_raise(d, v, 0, K)

        c = dmp_mul(d, dmp_sub(R, e, v, K), v, K)
        r = dmp_add(r, c, v, K)

        r = dmp_ground_trunc(r, p, v, K)

        D = dup_mul(D, [K.one, -a], K)
        D = dup_trunc(D, p, K)

    return r

def _collins_crt(r, R, P, p, K):
    """Wrapper of CRT for Collins's resultant algorithm. """
    return gf_int(gf_crt([r, R], [P, p], K), P*p)

@cythonized("u,v,n,m")
def dmp_zz_collins_resultant(f, g, u, K):
    """Collins's modular resultant algorithm in `Z[X]`. """

    n = dmp_degree(f, u)
    m = dmp_degree(g, u)

    if n < 0 or m < 0:
        return dmp_zero(u-1)

    A = dmp_max_norm(f, u, K)
    B = dmp_max_norm(g, u, K)

    a = dmp_ground_LC(f, u, K)
    b = dmp_ground_LC(g, u, K)

    v = u - 1

    B = K(2)*K.factorial(n+m)*A**m*B**n
    r, p, P = dmp_zero(v), K.one, K.one

    while P <= B:
        p = K(nextprime(p))

        while not (a % p) or not (b % p):
            p = K(nextprime(p))

        F = dmp_ground_trunc(f, p, u, K)
        G = dmp_ground_trunc(g, p, u, K)

        try:
            R = dmp_zz_modular_resultant(F, G, p, u, K)
        except HomomorphismFailed:
            continue

        if K.is_one(P):
            r = R
        else:
            r = dmp_apply_pairs(r, R, _collins_crt, (P, p, K), v, K)

        P *= p

    return r

@cythonized("u,n,m")
def dmp_qq_collins_resultant(f, g, u, K0):
    """Collins's modular resultant algorithm in `Q[X]`. """
    n = dmp_degree(f, u)
    m = dmp_degree(g, u)

    if n < 0 or m < 0:
        return dmp_zero(u-1)

    K1 = K0.get_ring()

    cf, f = dmp_clear_denoms(f, u, K0, K1)
    cg, g = dmp_clear_denoms(g, u, K0, K1)

    f = dmp_convert(f, u, K0, K1)
    g = dmp_convert(g, u, K0, K1)

    r = dmp_zz_collins_resultant(f, g, u, K1)
    r = dmp_convert(r, u-1, K1, K0)

    c = K0.convert(cf**m * cg**n, K1)

    return dmp_exquo_ground(r, c, u-1, K0)

@cythonized("u")
def dmp_resultant(f, g, u, K):
    """Computes resultant of two polynomials in `K[X]`. """
    if not u:
        return dup_resultant(f, g, K)

    if K.has_Field:
        if K.is_QQ and query('USE_COLLINS_RESULTANT'):
            return dmp_qq_collins_resultant(f, g, u, K)
    else:
        if K.is_ZZ and query('USE_COLLINS_RESULTANT'):
            return dmp_zz_collins_resultant(f, g, u, K)

    return dmp_prs_resultant(f, g, u, K)[0]

@cythonized("d,s")
def dup_discriminant(f, K):
    """Computes discriminant of a polynomial in `K[x]`. """
    d = dup_degree(f)

    if d <= 0:
        return K.zero
    else:
        s = (-1)**((d*(d-1)) // 2)
        c = dup_LC(f, K)

        r = dup_resultant(f, dup_diff(f, 1, K), K)

        return K.quo(r, c*K(s))

@cythonized("u,v,d,s")
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
        c = dmp_mul_ground(c, K(s), v, K)

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

@cythonized("u")
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
    elif query('USE_SIMPLIFY_GCD'):
        return _dmp_simplify_gcd(f, g, u, K)
    else:
        return None

@cythonized("u")
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
    elif query('USE_SIMPLIFY_GCD'):
        return _dmp_simplify_gcd(f, g, u, K)
    else:
        return None

@cythonized("u,v,df,dg")
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

@cythonized("u")
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

@cythonized("u")
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

def _dup_zz_gcd_interpolate(h, x, K):
    """Interpolate polynomial GCD from integer GCD. """
    f = []

    while h:
        g = h % x

        if g > x // 2:
            g -= x

        f.insert(0, g)
        h = (h-g) // x

    return f

@cythonized("i,df,dg")
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

            h = _dup_zz_gcd_interpolate(h, x, K)
            h = dup_primitive(h, K)[1]

            cff_, r = dup_div(f, h, K)

            if not r:
                cfg_, r = dup_div(g, h, K)

                if not r:
                    h = dup_mul_ground(h, gcd, K)
                    return h, cff_, cfg_

            cff = _dup_zz_gcd_interpolate(cff, x, K)

            h, r = dup_div(f, cff, K)

            if not r:
                cfg_, r = dup_div(g, h, K)

                if not r:
                    h = dup_mul_ground(h, gcd, K)
                    return h, cff, cfg_

            cfg = _dup_zz_gcd_interpolate(cfg, x, K)

            h, r = dup_div(g, cfg, K)

            if not r:
                cff_, r = dup_div(f, h, K)

                if not r:
                    h = dup_mul_ground(h, gcd, K)
                    return h, cff, cfg

        x = 73794*x * K.sqrt(K.sqrt(x)) // 27011

    raise HeuristicGCDFailed('no luck')

@cythonized("v")
def _dmp_zz_gcd_interpolate(h, x, v, K):
    """Interpolate polynomial GCD from integer GCD. """
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

@cythonized("u,v,i,dg,df")
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

            h = _dmp_zz_gcd_interpolate(h, x, v, K)
            h = dmp_ground_primitive(h, u, K)[1]

            cff_, r = dmp_div(f, h, u, K)

            if dmp_zero_p(r, u):
                cfg_, r = dmp_div(g, h, u, K)

                if dmp_zero_p(r, u):
                    h = dmp_mul_ground(h, gcd, u, K)
                    return h, cff_, cfg_

            cff = _dmp_zz_gcd_interpolate(cff, x, v, K)

            h, r = dmp_div(f, cff, u, K)

            if dmp_zero_p(r, u):
                cfg_, r = dmp_div(g, h, u, K)

                if dmp_zero_p(r, u):
                    h = dmp_mul_ground(h, gcd, u, K)
                    return h, cff, cfg_

            cfg = _dmp_zz_gcd_interpolate(cfg, x, v, K)

            h, r = dmp_div(g, cfg, u, K)

            if dmp_zero_p(r, u):
                cff_, r = dmp_div(f, h, u, K)

                if dmp_zero_p(r, u):
                    h = dmp_mul_ground(h, gcd, u, K)
                    return h, cff_, cfg

        x = 73794*x * K.sqrt(K.sqrt(x)) // 27011

    raise HeuristicGCDFailed('no luck')

def dup_qq_heu_gcd(f, g, K0):
    """Heuristic polynomial GCD in `Q[x]`. """
    result = _dup_ff_trivial_gcd(f, g, K0)

    if result is not None:
        return result

    K1 = K0.get_ring()

    cf, f = dup_clear_denoms(f, K0, K1)
    cg, g = dup_clear_denoms(g, K0, K1)

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

@cythonized("u")
def dmp_qq_heu_gcd(f, g, u, K0):
    """Heuristic polynomial GCD in `Q[X]`. """
    result = _dmp_ff_trivial_gcd(f, g, u, K0)

    if result is not None:
        return result

    K1 = K0.get_ring()

    cf, f = dmp_clear_denoms(f, u, K0, K1)
    cg, g = dmp_clear_denoms(g, u, K0, K1)

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

def dup_inner_gcd(f, g, K):
    """Computes polynomial GCD and cofactors of `f` and `g` in `K[x]`. """
    if K.has_Field or not K.is_Exact:
        if K.is_QQ and query('USE_HEU_GCD'):
            try:
                return dup_qq_heu_gcd(f, g, K)
            except HeuristicGCDFailed:
                pass

        return dup_ff_prs_gcd(f, g, K)
    else:
        if K.is_ZZ and query('USE_HEU_GCD'):
            try:
                return dup_zz_heu_gcd(f, g, K)
            except HeuristicGCDFailed:
                pass

        return dup_rr_prs_gcd(f, g, K)

@cythonized("u")
def _dmp_inner_gcd(f, g, u, K):
    """Helper function for `dmp_inner_gcd()`. """
    if K.has_Field or not K.is_Exact:
        if K.is_QQ and query('USE_HEU_GCD'):
            try:
                return dmp_qq_heu_gcd(f, g, u, K)
            except HeuristicGCDFailed:
                pass

        return dmp_ff_prs_gcd(f, g, u, K)
    else:
        if K.is_ZZ and query('USE_HEU_GCD'):
            try:
                 return dmp_zz_heu_gcd(f, g, u, K)
            except HeuristicGCDFailed:
                pass

        return dmp_rr_prs_gcd(f, g, u, K)

@cythonized("u")
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

@cythonized("u")
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

    return dup_monic(h, K)

def dup_lcm(f, g, K):
    """Computes polynomial LCM of `f` and `g` in `K[x]`. """
    if K.has_Field or not K.is_Exact:
        return dup_ff_lcm(f, g, K)
    else:
        return dup_rr_lcm(f, g, K)

@cythonized("u")
def dmp_rr_lcm(f, g, u, K):
    """Computes polynomial LCM over a ring in `K[X]`. """
    fc, f = dmp_ground_primitive(f, u, K)
    gc, g = dmp_ground_primitive(g, u, K)

    c = K.lcm(fc, gc)

    h = dmp_exquo(dmp_mul(f, g, u, K),
                  dmp_gcd(f, g, u, K), u, K)

    return dmp_mul_ground(h, c, u, K)

@cythonized("u")
def dmp_ff_lcm(f, g, u, K):
    """Computes polynomial LCM over a field in `K[X]`. """
    h = dmp_exquo(dmp_mul(f, g, u, K),
                  dmp_gcd(f, g, u, K), u, K)

    return dmp_ground_monic(h, u, K)

@cythonized("u")
def dmp_lcm(f, g, u, K):
    """Computes polynomial LCM of `f` and `g` in `K[X]`. """
    if not u:
        return dup_lcm(f, g, K)

    if K.has_Field or not K.is_Exact:
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

@cythonized("u")
def dmp_trunc(f, p, u, K):
    """Reduce `K[X]` polynomial modulo a polynomial `p` in `K[Y]`. """
    return dmp_strip([ dmp_rem(c, p, u-1, K) for c in f ], u)

@cythonized("u,v")
def dmp_ground_trunc(f, p, u, K):
    """Reduce `K[X]` polynomial modulo a constant `p` in `K`. """
    if not u:
        return dup_trunc(f, p, K)

    v = u-1

    return dmp_strip([ dmp_ground_trunc(c, p, v, K) for c in f ], u)

def dup_monic(f, K):
    """Divides all coefficients by `LC(f)` in `K[x]`. """
    if not f:
        return f

    lc = dup_LC(f, K)

    if K.is_one(lc):
        return f
    else:
        return dup_quo_ground(f, lc, K)

@cythonized("u")
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
    if K.has_Field or not K.is_Exact:
        return dup_ff_content(f, K)
    else:
        return dup_rr_content(f, K)

@cythonized("u,v")
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

@cythonized("u,v")
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

@cythonized("u")
def dmp_ff_ground_content(f, u, K):
    """Returns GCD of coefficients over a field. """
    if dmp_zero_p(f, u):
        return K.zero
    else:
        return K.one

@cythonized("u")
def dmp_ground_content(f, u, K):
    """Returns GCD of coefficients in `K[X]`. """
    if not u:
        return dup_content(f, K)

    if K.has_Field or not K.is_Exact:
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
    if K.has_Field or not K.is_Exact:
        return dup_ff_primitive(f, K)
    else:
        return dup_rr_primitive(f, K)

@cythonized("u,v")
def dmp_primitive(f, u, K):
    """Returns multivariate content and a primitive polynomial. """
    cont, v = dmp_content(f, u, K), u-1

    if dmp_zero_p(f, u) or dmp_one_p(cont, v, K):
        return cont, f
    else:
        return cont, [ dmp_exquo(c, cont, v, K) for c in f ]

@cythonized("u")
def dmp_rr_ground_primitive(f, u, K):
    """Returns content and a primitive polynomial over a ring. """
    cont = dmp_ground_content(f, u, K)

    if K.is_one(cont):
        return cont, f
    else:
        return cont, dmp_exquo_ground(f, cont, u, K)

@cythonized("u")
def dmp_ff_ground_primitive(f, u, K):
    """Returns content and a primitive polynomial over a ring. """
    if dmp_zero_p(f, u):
        return K.zero, f
    else:
        return K.one, f

@cythonized("u")
def dmp_ground_primitive(f, u, K):
    """Returns content and a primitive polynomial in `K[x]`. """
    if not u:
        return dup_primitive(f, K)

    if dmp_zero_p(f, u):
        return K.zero, f

    if K.has_Field or not K.is_Exact:
        return dmp_ff_ground_primitive(f, u, K)
    else:
        return dmp_rr_ground_primitive(f, u, K)

def dup_sqf_p(f, K):
    """Returns `True` if `f` is a square-free polynomial in `K[x]`. """
    if not f:
        return True
    else:
        return not dup_degree(dup_gcd(f, dup_diff(f, 1, K), K))

@cythonized("u")
def dmp_sqf_p(f, u, K):
    """Returns `True` if `f` is a square-free polynomial in `K[X]`. """
    if dmp_zero_p(f, u):
        return True
    else:
        return not dmp_degree(dmp_gcd(f, dmp_diff(f, 1, u, K), u, K), u)

@cythonized("s")
def dup_sqf_norm(f, K):
    """Square-free norm of `f` in `K[x]`, useful over algebraic domains. """
    if not K.is_Algebraic:
        raise DomainError("ground domain must be algebraic")

    s, g = 0, dmp_raise(K.mod.rep, 1, 0, K.dom)

    while True:
        h, _ = dmp_inject(f, 0, K, front=True)
        r = dmp_resultant(g, h, 1, K.dom)

        if dup_sqf_p(r, K.dom):
            break
        else:
            f, s = dup_taylor(f, -K.unit, K), s+1

    return s, f, r

@cythonized("s,u")
def dmp_sqf_norm(f, u, K):
    """Square-free norm of `f` in `K[X]`, useful over algebraic domains. """
    if not u:
        return dup_sqf_norm(f, K)

    if not K.is_Algebraic:
        raise DomainError("ground domain must be algebraic")

    g = dmp_raise(K.mod.rep, u+1, 0, K.dom)
    F = dmp_raise([K.one,-K.unit], u, 0, K)

    s = 0

    while True:
        h, _ = dmp_inject(f, u, K, front=True)
        r = dmp_resultant(g, h, u+1, K.dom)

        if dmp_sqf_p(r, u, K.dom):
            break
        else:
            f, s = dmp_compose(f, F, u, K), s+1

    return s, f, r

def dup_sqf_part(f, K):
    """Returns square-free part of a polynomial in `K[x]`. """
    if not f:
        return f

    if K.is_negative(dup_LC(f, K)):
        f = dup_neg(f, K)

    gcd = dup_gcd(f, dup_diff(f, 1, K), K)
    sqf = dup_exquo(f, gcd, K)

    if K.has_Field or not K.is_Exact:
        return dup_monic(sqf, K)
    else:
        return dup_primitive(sqf, K)[1]

@cythonized("u")
def dmp_sqf_part(f, u, K):
    """Returns square-free part of a polynomial in `K[X]`. """
    if dmp_zero_p(f, u):
        return f

    if K.is_negative(dmp_ground_LC(f, u, K)):
        f = dmp_neg(f, u, K)

    gcd = dmp_gcd(f, dmp_diff(f, 1, u, K), u, K)
    sqf = dmp_exquo(f, gcd, u, K)

    if K.has_Field or not K.is_Exact:
        return dmp_ground_monic(sqf, u, K)
    else:
        return dmp_ground_primitive(sqf, u, K)[1]

@cythonized("i")
def dup_sqf_list(f, K, all=False):
    """Returns square-free decomposition of a polynomial in `K[x]`. """
    if K.has_Field or not K.is_Exact:
        coeff = dup_LC(f, K)
        f = dup_monic(f, K)
    else:
        coeff, f = dup_primitive(f, K)

        if K.is_negative(dup_LC(f, K)):
            f = dup_neg(f, K)
            coeff = -coeff

    if dup_degree(f) <= 0:
        return coeff, []

    result, i = [], 1

    h = dup_diff(f, 1, K)
    g, p, q = dup_inner_gcd(f, h, K)

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

    return coeff, result

def dup_sqf_list_include(f, K, all=False):
    """Returns square-free decomposition of a polynomial in `K[x]`. """
    coeff, factors = dup_sqf_list(f, K, all)

    if not factors:
        return [(dup_strip([coeff]), 1)]
    else:
        g = dup_mul_ground(factors[0][0], coeff, K)
        return [(g, factors[0][1])] + factors[1:]

@cythonized("u,i")
def dmp_sqf_list(f, u, K, all=False):
    """Returns square-free decomposition of a polynomial in `K[X]`. """
    if not u:
        return dup_sqf_list(f, K, all)

    if K.has_Field or not K.is_Exact:
        coeff = dmp_ground_LC(f, u, K)
        f = dmp_ground_monic(f, u, K)
    else:
        coeff, f = dmp_ground_primitive(f, u, K)

        if K.is_negative(dmp_ground_LC(f, u, K)):
            f = dmp_neg(f, u, K)
            coeff = -coeff

    if dmp_degree(f, u) <= 0:
        return coeff, []

    result, i = [], 1

    h = dmp_diff(f, 1, u, K)
    g, p, q = dmp_inner_gcd(f, h, u, K)

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

    return coeff, result

@cythonized("u")
def dmp_sqf_list_include(f, u, K, all=False):
    """Returns square-free decomposition of a polynomial in `K[x]`. """
    if not u:
        return dup_sqf_list_include(f, K, all)

    coeff, factors = dmp_sqf_list(f, u, K, all)

    if not factors:
        return [(dmp_ground(coeff, u), 1)]
    else:
        g = dmp_mul_ground(factors[0][0], coeff, u, K)
        return [(g, factors[0][1])] + factors[1:]

def dup_extract(f, g, K):
    """Extracts common content from a pair of polynomials in `K[x]`. """
    fc = dup_content(f, K)
    gc = dup_content(g, K)

    gcd = K.gcd(fc, gc)

    if not K.is_one(gcd):
        f = dup_exquo_ground(f, gcd, K)
        g = dup_exquo_ground(g, gcd, K)

    return gcd, f, g

@cythonized("u")
def dmp_ground_extract(f, g, u, K):
    """Extracts common content from a pair of polynomials in `K[X]`. """
    fc = dmp_ground_content(f, u, K)
    gc = dmp_ground_content(g, u, K)

    gcd = K.gcd(fc, gc)

    if not K.is_one(gcd):
        f = dmp_exquo_ground(f, gcd, u, K)
        g = dmp_exquo_ground(g, gcd, u, K)

    return gcd, f, g

def dup_real_imag(f, K):
    """Return bivariate polynomials f1 and f2, such that f = f1 + f2*I. """
    if not K.is_ZZ and not K.is_QQ:
        raise DomainError("computing real and imaginary parts is not supported over %s" % K)

    f1 = dmp_zero(1)
    f2 = dmp_zero(1)

    if not f:
        return f1, f2

    g = [[[K.one, K.zero]], [[K.one], []]]
    h = dmp_ground(f[0], 2)

    for c in f[1:]:
        h = dmp_mul(h, g, 2, K)
        h = dmp_add_term(h, dmp_ground(c, 1), 0, 2, K)

    H = dup_to_raw_dict(h)

    for k, h in H.iteritems():
        m = k % 4

        if not m:
            f1 = dmp_add(f1, h, 1, K)
        elif m == 1:
            f2 = dmp_add(f2, h, 1, K)
        elif m == 2:
            f1 = dmp_sub(f1, h, 1, K)
        else:
            f2 = dmp_sub(f2, h, 1, K)

    return f1, f2

def dup_mirror(f, K):
    """Evaluate efficiently composition `f(-x)` in `K[x]`. """
    f, n, a = list(f), dup_degree(f), -K.one

    for i in xrange(n-1, -1, -1):
        f[i], a = a*f[i], -a

    return f

def dup_scale(f, a, K):
    """Evaluate efficiently composition `f(a*x)` in `K[x]`. """
    f, n, b = list(f), dup_degree(f), a

    for i in xrange(n-1, -1, -1):
        f[i], b = b*f[i], b*a

    return f

def dup_taylor(f, a, K):
    """Evaluate efficiently Taylor shift `f(x + a)` in `K[x]`. """
    f, n = list(f), dup_degree(f)

    for i in xrange(n, 0, -1):
        for j in xrange(0, i):
            f[j+1] += a*f[j]

    return f

def dup_transform(f, p, q, K):
    """Evaluate functional transformation `q**n * f(p/q)` in `K[x]`. """
    if not f:
        return []

    h, Q = [f[0]], [[K.one]]

    for i in xrange(0, dup_degree(f)):
        Q.append(dup_mul(Q[-1], q, K))

    for c, q in zip(f[1:], Q[1:]):
        h = dup_mul(h, p, K)
        q = dup_mul_ground(q, c, K)
        h = dup_add(h, q, K)

    return h

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

@cythonized("u")
def dmp_compose(f, g, u, K):
    """Evaluate functional composition `f(g)` in `K[X]`. """
    if not u:
        return dup_compose(f, g, K)

    if dmp_zero_p(f, u):
        return f

    h = [f[0]]

    for c in f[1:]:
        h = dmp_mul(h, g, u, K)
        h = dmp_add_term(h, c, 0, u, K)

    return h

@cythonized("s,n,r,i,j")
def _dup_right_decompose(f, s, K):
    """XXX"""
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

        g[s-i] = K.exquo(coeff, i*r*lc)

    return dup_from_raw_dict(g, K)

@cythonized("i")
def _dup_left_decompose(f, h, K):
    """XXX"""
    g, i = {}, 0

    while f:
        q, r = dup_div(f, h, K)

        if dup_degree(r) > 0:
            return None
        else:
            g[i] = dup_LC(r, K)
            f, i = q, i + 1

    return dup_from_raw_dict(g, K)


@cythonized("df,s")
def _dup_decompose(f, K):
    """XXX"""
    df = dup_degree(f)

    for s in xrange(2, df):
        if df % s != 0:
            continue

        h = _dup_right_decompose(f, s, K)

        if h is not None:
            g = _dup_left_decompose(f, h, K)

            if g is not None:
                return g, h

    return None

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
    F = []

    while True:
        result = _dup_decompose(f, K)

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
    if not (K.has_Field or not K.is_Exact):
        raise DomainError("can't compute Sturm sequence over %s" % K)

    f = dup_sqf_part(f, K)

    sturm = [f, dup_diff(f, 1, K)]

    while sturm[-1]:
        s = dup_rem(sturm[-2], sturm[-1], K)
        sturm.append(dup_neg(s, K))

    return sturm[:-1]

@cythonized("u")
def dmp_lift(f, u, K):
    """Convert algebraic coefficients to integers in `K[X]`. """
    if not K.is_Algebraic:
        raise DomainError('computation can be done only in an algebraic domain')

    F, monoms, polys = dmp_to_dict(f, u), [], []

    for monom, coeff in F.iteritems():
        if not coeff.is_ground:
            monoms.append(monom)

    perms = variations([-1, 1], len(monoms), repetition=True)

    for perm in perms:
        G = dict(F)

        for sign, monom in zip(perm, monoms):
            if sign == -1:
                G[monom] = -G[monom]

        polys.append(dmp_from_dict(G, u, K))

    return dmp_convert(dmp_expand(polys, u, K), u, K, K.dom)

def dup_sign_variations(f, K):
    """Compute the number of sign variations of `f` in `K[x]`. """
    prev, k = K.zero, 0

    for coeff in f:
        if coeff*prev < 0:
            k += 1

        if coeff:
            prev = coeff

    return k

