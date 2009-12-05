"""
This module implements computation of elementary transcendental
functions (powers, logarithms, trigonometric and hyperbolic
functions, inverse trigonometric and hyperbolic) for real
floating-point numbers.

For complex and interval implementations of the same functions,
see libmpc and libmpi.

"""

import math
from bisect import bisect

from settings import (
    MP_BASE, MP_ZERO, MP_ONE, MP_TWO, MP_FIVE, MODE,
    round_floor, round_ceiling, round_down, round_up,
    round_nearest, round_fast,
)

from libmpf import (
    ComplexResult,
    bitcount, bctable, lshift, rshift, giant_steps, sqrt_fixed,
    from_int, to_int, from_man_exp, to_fixed, to_float, from_float,
    normalize,
    fzero, fone, fnone, fhalf, finf, fninf, fnan,
    mpf_cmp, mpf_sign, mpf_abs,
    mpf_pos, mpf_neg, mpf_add, mpf_sub, mpf_mul, mpf_div, mpf_shift,
    mpf_rdiv_int, mpf_pow_int, mpf_sqrt,
    reciprocal_rnd, negative_rnd, mpf_perturb,
    isqrt_fast
)

from libintmath import ifib


#----------------------------------------------------------------------------#
#                                                                            #
#                   Elementary mathematical constants                        #
#                                                                            #
#----------------------------------------------------------------------------#

def constant_memo(f):
    """
    Decorator for caching computed values of mathematical
    constants. This decorator should be applied to a
    function taking a single argument prec as input and
    returning a fixed-point value with the given precision.
    """
    f.memo_prec = -1
    f.memo_val = None
    def g(prec, **kwargs):
        memo_prec = f.memo_prec
        if prec <= memo_prec:
            return f.memo_val >> (memo_prec-prec)
        newprec = int(prec*1.05+10)
        f.memo_val = f(newprec, **kwargs)
        f.memo_prec = newprec
        return f.memo_val >> (newprec-prec)
    g.__name__ = f.__name__
    g.__doc__ = f.__doc__
    return g

def def_mpf_constant(fixed):
    """
    Create a function that computes the mpf value for a mathematical
    constant, given a function that computes the fixed-point value.

    Assumptions: the constant is positive and has magnitude ~= 1;
    the fixed-point function rounds to floor.
    """
    def f(prec, rnd=round_fast):
        wp = prec + 20
        v = fixed(wp)
        if rnd in (round_up, round_ceiling):
            v += 1
        return normalize(0, v, -wp, bitcount(v), prec, rnd)
    f.__doc__ = fixed.__doc__
    return f

def bsp_acot(q, a, b, hyperbolic):
    if b - a == 1:
        a1 = MP_BASE(2*a + 3)
        if hyperbolic or a&1:
            return MP_ONE, a1 * q**2, a1
        else:
            return -MP_ONE, a1 * q**2, a1
    m = (a+b)//2
    p1, q1, r1 = bsp_acot(q, a, m, hyperbolic)
    p2, q2, r2 = bsp_acot(q, m, b, hyperbolic)
    return q2*p1 + r1*p2, q1*q2, r1*r2

# the acoth(x) series converges like the geometric series for x^2
# N = ceil(p*log(2)/(2*log(x)))
def acot_fixed(a, prec, hyperbolic):
    """
    Compute acot(a) or acoth(a) for an integer a with binary splitting; see
    http://numbers.computation.free.fr/Constants/Algorithms/splitting.html
    """
    N = int(0.35 * prec/math.log(a) + 20)
    p, q, r = bsp_acot(a, 0,N, hyperbolic)
    return ((p+q)<<prec)//(q*a)

def machin(coefs, prec, hyperbolic=False):
    """
    Evaluate a Machin-like formula, i.e., a linear combination of
    acot(n) or acoth(n) for specific integer values of n, using fixed-
    point arithmetic. The input should be a list [(c, n), ...], giving
    c*acot[h](n) + ...
    """
    extraprec = 10
    s = MP_ZERO
    for a, b in coefs:
        s += MP_BASE(a) * acot_fixed(MP_BASE(b), prec+extraprec, hyperbolic)
    return (s >> extraprec)

# Logarithms of integers are needed for various computations involving
# logarithms, powers, radix conversion, etc

@constant_memo
def ln2_fixed(prec):
    """
    Computes ln(2). This is done with a hyperbolic Machin-type formula,
    with binary splitting at high precision.
    """
    return machin([(18, 26), (-2, 4801), (8, 8749)], prec, True)

@constant_memo
def ln10_fixed(prec):
    """
    Computes ln(10). This is done with a hyperbolic Machin-type formula.
    """
    return machin([(46, 31), (34, 49), (20, 161)], prec, True)


"""
For computation of pi, we use the Chudnovsky series:

             oo
             ___        k
      1     \       (-1)  (6 k)! (A + B k)
    ----- =  )     -----------------------
    12 pi   /___               3  3k+3/2
                    (3 k)! (k!)  C
            k = 0

where A, B, and C are certain integer constants. This series adds roughly
14 digits per term. Note that C^(3/2) can be extracted so that the
series contains only rational terms. This makes binary splitting very
efficient.

The recurrence formulas for the binary splitting were taken from
ftp://ftp.gmplib.org/pub/src/gmp-chudnovsky.c

Previously, Machin's formula was used at low precision and the AGM iteration
was used at high precision. However, the Chudnovsky series is essentially as
fast as the Machin formula at low precision and in practice about 3x faster
than the AGM at high precision (despite theoretically having a worse
asymptotic complexity), so there is no reason not to use it in all cases.

"""

# Constants in Chudnovsky's series
CHUD_A = MP_BASE(13591409)
CHUD_B = MP_BASE(545140134)
CHUD_C = MP_BASE(640320)
CHUD_D = MP_BASE(12)

def bs_chudnovsky(a, b, level, verbose):
    """
    Computes the sum from a to b of the series in the Chudnovsky
    formula. Returns g, p, q where p/q is the sum as an exact
    fraction and g is a temporary value used to save work
    for recursive calls.
    """
    if b-a == 1:
        g = MP_BASE((6*b-5)*(2*b-1)*(6*b-1))
        p = b**3 * CHUD_C**3 // 24
        q = (-1)**b * g * (CHUD_A+CHUD_B*b)
    else:
        if verbose and level < 4:
            print "  binary splitting", a, b
        mid = (a+b)//2
        g1, p1, q1 = bs_chudnovsky(a, mid, level+1, verbose)
        g2, p2, q2 = bs_chudnovsky(mid, b, level+1, verbose)
        p = p1*p2
        g = g1*g2
        q = q1*p2 + q2*g1
    return g, p, q

@constant_memo
def pi_fixed(prec, verbose=False, verbose_base=None):
    """
    Compute floor(pi * 2**prec) as a big integer.

    This is done using Chudnovsky's series (see comments in
    libelefun.py for details).
    """
    # The Chudnovsky series gives 14.18 digits per term
    N = int(prec/3.3219280948/14.181647462 + 2)
    if verbose:
        print "binary splitting with N =", N
    g, p, q = bs_chudnovsky(0, N, 0, verbose)
    sqrtC = isqrt_fast(CHUD_C<<(2*prec))
    v = p*CHUD_C*sqrtC//((q+CHUD_A*p)*CHUD_D)
    return v

def degree_fixed(prec):
    return pi_fixed(prec)//180

def bspe(a, b):
    """
    Sum series for exp(1)-1 between a, b, returning the result
    as an exact fraction (p, q).
    """
    if b-a == 1:
        return MP_ONE, MP_BASE(b)
    m = (a+b)//2
    p1, q1 = bspe(a, m)
    p2, q2 = bspe(m, b)
    return p1*q2+p2, q1*q2

@constant_memo
def e_fixed(prec):
    """
    Computes exp(1). This is done using the ordinary Taylor series for
    exp, with binary splitting. For a description of the algorithm,
    see:

        http://numbers.computation.free.fr/Constants/
            Algorithms/splitting.html
    """
    # Slight overestimate of N needed for 1/N! < 2**(-prec)
    # This could be tightened for large N.
    N = int(1.1*prec/math.log(prec) + 20)
    p, q = bspe(0,N)
    return ((p+q)<<prec)//q

@constant_memo
def phi_fixed(prec):
    """
    Computes the golden ratio, (1+sqrt(5))/2
    """
    prec += 10
    a = isqrt_fast(MP_FIVE<<(2*prec)) + (MP_ONE << prec)
    return a >> 11

mpf_phi    = def_mpf_constant(phi_fixed)
mpf_pi     = def_mpf_constant(pi_fixed)
mpf_e      = def_mpf_constant(e_fixed)
mpf_degree = def_mpf_constant(degree_fixed)
mpf_ln2    = def_mpf_constant(ln2_fixed)
mpf_ln10   = def_mpf_constant(ln10_fixed)


#----------------------------------------------------------------------------#
#                                                                            #
#                                    Powers                                  #
#                                                                            #
#----------------------------------------------------------------------------#

def mpf_pow(s, t, prec, rnd=round_fast):
    """
    Compute s**t. Raises ComplexResult if s is negative and t is
    fractional.
    """
    ssign, sman, sexp, sbc = s
    tsign, tman, texp, tbc = t
    if ssign and texp < 0:
        raise ComplexResult("negative number raised to a fractional power")
    if texp >= 0:
        return mpf_pow_int(s, (-1)**tsign * (tman<<texp), prec, rnd)
    # s**(n/2) = sqrt(s)**n
    if texp == -1:
        if tman == 1:
            if tsign:
                return mpf_div(fone, mpf_sqrt(s, prec+10,
                    reciprocal_rnd[rnd]), prec, rnd)
            return mpf_sqrt(s, prec, rnd)
        else:
            if tsign:
                return mpf_pow_int(mpf_sqrt(s, prec+10,
                    reciprocal_rnd[rnd]), -tman, prec, rnd)
            return mpf_pow_int(mpf_sqrt(s, prec+10, rnd), tman, prec, rnd)
    # General formula: s**t = exp(t*log(s))
    # TODO: handle rnd direction of the logarithm carefully
    c = mpf_log(s, prec+10, rnd)
    return mpf_exp(mpf_mul(t, c), prec, rnd)

def int_pow_fixed(y, n, prec):
    """n-th power of a fixed point number with precision prec

       Returns the power in the form man, exp,
       man * 2**exp ~= y**n
    """
    if n == 2:
        return (y*y), 0
    bc = bitcount(y)
    exp = 0
    workprec = 2 * (prec + 4*bitcount(n) + 4)
    _, pm, pe, pbc = fone
    while 1:
        if n & 1:
            pm = pm*y
            pe = pe+exp
            pbc += bc - 2
            pbc = pbc + bctable[int(pm >> pbc)]
            if pbc > workprec:
                pm = pm >> (pbc-workprec)
                pe += pbc - workprec
                pbc = workprec
            n -= 1
            if not n:
                break
        y = y*y
        exp = exp+exp
        bc = bc + bc - 2
        bc = bc + bctable[int(y >> bc)]
        if bc > workprec:
            y = y >> (bc-workprec)
            exp += bc - workprec
            bc = workprec
        n = n // 2
    return pm, pe

# froot(s, n, prec, rnd) computes the real n-th root of a
# positive mpf tuple s.
# To compute the root we start from a 50-bit estimate for r
# generated with ordinary floating-point arithmetic, and then refine
# the value to full accuracy using the iteration

#            1  /                     y       \
#   r     = --- | (n-1)  * r   +  ----------  |
#    n+1     n  \           n     r_n**(n-1)  /

# which is simply Newton's method applied to the equation r**n = y.
# With giant_steps(start, prec+extra) = [p0,...,pm, prec+extra]
# and y = man * 2**-shift  one has
# (man * 2**exp)**(1/n) =
# y**(1/n) * 2**(start-prec/n) * 2**(p0-start) * ... * 2**(prec+extra-pm) *
# 2**((exp+shift-(n-1)*prec)/n -extra))
# The last factor is accounted for in the last line of froot.

def nthroot_fixed(y, n, prec, exp1):
    start = 50
    try:
        y1 = rshift(y, prec - n*start)
        r = MP_BASE(int(y1**(1.0/n)))
    except OverflowError:
        y1 = from_int(y1, start)
        fn = from_int(n)
        fn = mpf_rdiv_int(1, fn, start)
        r = mpf_pow(y1, fn, start)
        r = to_int(r)
    extra = 10
    extra1 = n
    prevp = start
    for p in giant_steps(start, prec+extra):
        pm, pe = int_pow_fixed(r, n-1, prevp)
        r2 = rshift(pm, (n-1)*prevp - p - pe - extra1)
        B = lshift(y, 2*p-prec+extra1)//r2
        r = (B + (n-1) * lshift(r, p-prevp))//n
        prevp = p
    return r

def mpf_nthroot(s, n, prec, rnd=round_fast):
    """nth-root of a positive number

    Use the Newton method when faster, otherwise use x**(1/n)
    """
    sign, man, exp, bc = s
    if sign:
        raise ComplexResult("nth root of a negative number")
    if not man:
        if s == fnan:
            return fnan
        if s == fzero:
            if n > 0:
                return fzero
            if n == 0:
                return fone
            return finf
        # Infinity
        if not n:
            return fnan
        if n < 0:
            return fzero
        return finf
    flag_inverse = False
    if n < 2:
        if n == 0:
            return fone
        if n == 1:
            return mpf_pos(s, prec, rnd)
        if n == -1:
            return mpf_div(fone, s, prec, rnd)
        # n < 0
        rnd = reciprocal_rnd[rnd]
        flag_inverse = True
        extra_inverse = 5
        prec += extra_inverse
        n = -n
    if n > 20 and (n >= 20000 or prec < int(233 + 28.3 * n**0.62)):
        prec2 = prec + 10
        fn = from_int(n)
        nth = mpf_rdiv_int(1, fn, prec2)
        r = mpf_pow(s, nth, prec2, rnd)
        s = normalize(r[0], r[1], r[2], r[3], prec, rnd)
        if flag_inverse:
            return mpf_div(fone, s, prec-extra_inverse, rnd)
        else:
            return s
    # Convert to a fixed-point number with prec2 bits.
    prec2 = prec + 2*n - (prec%n)
    # a few tests indicate that
    # for 10 < n < 10**4 a bit more precision is needed
    if n > 10:
        prec2 += prec2//10
        prec2 = prec2 - prec2%n
    # Mantissa may have more bits than we need. Trim it down.
    shift = bc - prec2
    # Adjust exponents to make prec2 and exp+shift multiples of n.
    sign1 = 0
    es = exp+shift
    if es < 0:
        sign1 = 1
        es = -es
    if sign1:
        shift += es%n
    else:
        shift -= es%n
    man = rshift(man, shift)
    extra = 10
    exp1 = ((exp+shift-(n-1)*prec2)//n) - extra
    rnd_shift = 0
    if flag_inverse:
        if rnd == 'u' or rnd == 'c':
            rnd_shift = 1
    else:
        if rnd == 'd' or rnd == 'f':
            rnd_shift = 1
    man = nthroot_fixed(man+rnd_shift, n, prec2, exp1)
    s = from_man_exp(man, exp1, prec, rnd)
    if flag_inverse:
        return mpf_div(fone, s, prec-extra_inverse, rnd)
    else:
        return s

def mpf_cbrt(s, prec, rnd=round_fast):
    """cubic root of a positive number"""
    return mpf_nthroot(s, 3, prec, rnd)

#----------------------------------------------------------------------------#
#                                                                            #
#                                Logarithms                                  #
#                                                                            #
#----------------------------------------------------------------------------#

# Fast sequential integer logarithms are required for various series
# computations related to zeta functions, so we cache them
# TODO: can this be done better?
log_int_cache = {}

def log_int_fixed(n, prec):
    if n < 2:
        return MP_ZERO
    cache = log_int_cache.get(prec)
    if cache and (n in cache):
        return cache[n]
    if cache:
        L = cache[max(cache)]
    else:
        cache = log_int_cache[prec] = {}
        L = cache[2] = ln2_fixed(prec)
    one = MP_ONE << prec
    for p in xrange(max(cache)+1, n+1):
        s = 0
        u = one
        k = 1
        a = (2*p-1)**2
        while u:
            s += u // k
            u //= a
            k += 2
        L += 2*s//(2*p-1)
        cache[p] = L
    return cache[n]

# Use Taylor series with caching up to this prec
LOG_TAYLOR_PREC = 2500

# Cache log values in steps of size 2^-N
LOG_TAYLOR_SHIFT = 9

# prec/size ratio of x for fastest convergence in AGM formula
LOG_AGM_MAG_PREC_RATIO = 20

log_taylor_cache = {}

# ~= next power of two + 20
cache_prec_steps = [22,22]
for k in xrange(1, bitcount(LOG_TAYLOR_PREC)+1):
    cache_prec_steps += [min(2**k,LOG_TAYLOR_PREC)+20] * 2**(k-1)

def agm_fixed(a, b, prec):
    """
    Fixed-point computation of agm(a,b), assuming
    a, b both close to unit magnitude.
    """
    i = 0
    while 1:
        anew = (a+b)>>1
        if i > 4 and abs(a-anew) < 8:
            return a
        b = isqrt_fast(a*b)
        a = anew
        i += 1
    return a

def log_agm(x, prec):
    """
    Fixed-point computation of -log(x) = log(1/x), suitable
    for large precision. It is required that 0 < x < 1. The
    algorithm used is the Sasaki-Kanada formula

        -log(x) = pi/agm(theta2(x)^2,theta3(x)^2). [1]

    For faster convergence in the theta functions, x should
    be chosen closer to 0.

    Guard bits must be added by the caller.

    HYPOTHESIS: if x = 2^(-n), n bits need to be added to
    account for the truncation to a fixed-point number,
    and this is the only significant cancellation error.

    The number of bits lost to roundoff is small and can be
    considered constant.

    [1] Richard P. Brent, "Fast Algorithms for High-Precision
        Computation of Elementary Functions (extended abstract)",
        http://wwwmaths.anu.edu.au/~brent/pd/RNC7-Brent.pdf

    """
    x2 = (x*x) >> prec
    # Compute jtheta2(x)**2
    s = a = b = x2
    while a:
        b = (b*x2) >> prec
        a = (a*b) >> prec
        s += a
    s += (MP_ONE<<prec)
    s = (s*s)>>(prec-2)
    s = (s*isqrt_fast(x<<prec))>>prec
    # Compute jtheta3(x)**2
    t = a = b = x
    while a:
        b = (b*x2) >> prec
        a = (a*b) >> prec
        t += a
    t = (MP_ONE<<prec) + (t<<1)
    t = (t*t)>>prec
    # Final formula
    p = agm_fixed(s, t, prec)
    return (pi_fixed(prec) << prec) // p

def log_taylor(x, prec, r=0):
    """
    Fixed-point calculation of log(x). It is assumed that x is close
    enough to 1 for the Taylor series to converge quickly. Convergence
    can be improved by specifying r > 0 to compute
    log(x^(1/2^r))*2^r, at the cost of performing r square roots.

    The caller must provide sufficient guard bits.
    """
    for i in xrange(r):
        x = isqrt_fast(x<<prec)
    one = MP_ONE << prec
    v = ((x-one)<<prec)//(x+one)
    sign = v < 0
    if sign:
        v = -v
    v2 = (v*v) >> prec
    v4 = (v2*v2) >> prec
    s0 = v
    s1 = v//3
    v = (v*v4) >> prec
    k = 5
    while v:
        s0 += v // k
        k += 2
        s1 += v // k
        v = (v*v4) >> prec
        k += 2
    s1 = (s1*v2) >> prec
    s = (s0+s1) << (1+r)
    if sign:
        return -s
    return s

def log_taylor_cached(x, prec):
    """
    Fixed-point computation of log(x), assuming x in (0.5, 2)
    and prec <= LOG_TAYLOR_PREC.
    """
    n = x >> (prec-LOG_TAYLOR_SHIFT)
    cached_prec = cache_prec_steps[prec]
    dprec = cached_prec - prec
    if (n, cached_prec) in log_taylor_cache:
        a, log_a = log_taylor_cache[n, cached_prec]
    else:
        a = n << (cached_prec - LOG_TAYLOR_SHIFT)
        log_a = log_taylor(a, cached_prec, 8)
        log_taylor_cache[n, cached_prec] = (a, log_a)
    a >>= dprec
    log_a >>= dprec
    u = ((x - a) << prec) // a
    v = (u << prec) // ((MP_TWO << prec) + u)
    v2 = (v*v) >> prec
    v4 = (v2*v2) >> prec
    s0 = v
    s1 = v//3
    v = (v*v4) >> prec
    k = 5
    while v:
        s0 += v//k
        k += 2
        s1 += v//k
        v = (v*v4) >> prec
        k += 2
    s1 = (s1*v2) >> prec
    s = (s0+s1) << 1
    return log_a + s

def mpf_log(x, prec, rnd=round_fast):
    """
    Compute the natural logarithm of the mpf value x. If x is negative,
    ComplexResult is raised.
    """
    sign, man, exp, bc = x
    #------------------------------------------------------------------
    # Handle special values
    if not man:
        if x == fzero: return fninf
        if x == finf: return finf
        if x == fnan: return fnan
    if sign:
        raise ComplexResult("logarithm of a negative number")
    wp = prec + 20
    #------------------------------------------------------------------
    # Handle log(2^n) = log(n)*2.
    # Here we catch the only possible exact value, log(1) = 0
    if man == 1:
        if not exp:
            return fzero
        return from_man_exp(exp*ln2_fixed(wp), -wp, prec, rnd)
    mag = exp+bc
    abs_mag = abs(mag)
    #------------------------------------------------------------------
    # Handle x = 1+eps, where log(x) ~ x. We need to check for
    # cancellation when moving to fixed-point math and compensate
    # by increasing the precision. Note that abs_mag in (0, 1) <=>
    # 0.5 < x < 2 and x != 1
    if abs_mag <= 1:
        # Calculate t = x-1 to measure distance from 1 in bits
        tsign = 1-abs_mag
        if tsign:
            tman = (MP_ONE<<bc) - man
        else:
            tman = man - (MP_ONE<<(bc-1))
        tbc = bitcount(tman)
        cancellation = bc - tbc
        if cancellation > wp:
            t = normalize(tsign, tman, abs_mag-bc, tbc, tbc, 'n')
            return mpf_perturb(t, tsign, prec, rnd)
        else:
            wp += cancellation
        # TODO: if close enough to 1, we could use Taylor series
        # even in the AGM precision range, since the Taylor series
        # converges rapidly
    #------------------------------------------------------------------
    # Another special case:
    # n*log(2) is a good enough approximation
    if abs_mag > 10000:
        if bitcount(abs_mag) > wp:
            return from_man_exp(exp*ln2_fixed(wp), -wp, prec, rnd)
    #------------------------------------------------------------------
    # General case.
    # Perform argument reduction using log(x) = log(x*2^n) - n*log(2):
    # If we are in the Taylor precision range, choose magnitude 0 or 1.
    # If we are in the AGM precision range, choose magnitude -m for
    # some large m; benchmarking on one machine showed m = prec/20 to be
    # optimal between 1000 and 100,000 digits.
    if wp <= LOG_TAYLOR_PREC:
        m = log_taylor_cached(lshift(man, wp-bc), wp)
        if mag:
            m += mag*ln2_fixed(wp)
    else:
        optimal_mag = -wp//LOG_AGM_MAG_PREC_RATIO
        n = optimal_mag - mag
        x = mpf_shift(x, n)
        wp += (-optimal_mag)
        m = -log_agm(to_fixed(x, wp), wp)
        m -= n*ln2_fixed(wp)
    return from_man_exp(m, -wp, prec, rnd)

def mpf_log_hypot(a, b, prec, rnd):
    """
    Computes log(sqrt(a^2+b^2)) accurately.
    """
    # If either a or b is inf/nan/0, assume it to be a
    if not b[1]:
        a, b = b, a
    # a is inf/nan/0
    if not a[1]:
        # both are inf/nan/0
        if not b[1]:
            if a == b == fzero:
                return fninf
            if fnan in (a, b):
                return fnan
            # at least one term is (+/- inf)^2
            return finf
        # only a is inf/nan/0
        if a == fzero:
            # log(sqrt(0+b^2)) = log(|b|)
            return mpf_log(mpf_abs(b), prec, rnd)
        if a == fnan:
            return fnan
        return finf
    # Exact
    a2 = mpf_mul(a,a)
    b2 = mpf_mul(b,b)
    extra = 20
    # Not exact
    h2 = mpf_add(a2, b2, prec+extra)
    cancelled = mpf_add(h2, fnone, 10)
    mag_cancelled = cancelled[2]+cancelled[3]
    # Just redo the sum exactly if necessary (could be smarter
    # and avoid memory allocation when a or b is precisely 1
    # and the other is tiny...)
    if cancelled == fzero or mag_cancelled < -extra//2:
        h2 = mpf_add(a2, b2, prec+extra-min(a2[2],b2[2]))
    return mpf_shift(mpf_log(h2, prec, rnd), -1)


#----------------------------------------------------------------------------#
#                                                                            #
#                           Exponential function                             #
#                                                                            #
#----------------------------------------------------------------------------#

# The exponential function has a rapidly convergent Maclaurin series:
#
#     exp(x) = 1 + x + x**2/2! + x**3/3! + x**4/4! + ...
#
# The series can be summed very easily using fixed-point arithmetic.
# The convergence can be improved further, using a trick due to
# Richard P. Brent: instead of computing exp(x) directly, we choose a
# small integer r (say, r=10) and compute exp(x/2**r)**(2**r).

# The optimal value for r depends on the Python platform, the magnitude
# of x and the target precision, and has to be estimated from
# experimental timings. One test with x ~= 0.3 showed that
# r = 2.2*prec**0.42 gave a good fit to the optimal values for r for
# prec between 1 and 10000 bits, on one particular machine.

# This optimization makes the summation about twice as fast at
# low precision levels and much faster at high precision
# (roughly five times faster at 1000 decimal digits).

# If |x| is very large, we first rewrite it as t + n*log(2) with the
# integer n chosen such that |t| <= log(2), and then calculate
# exp(x) as exp(t)*(2**n), using the Maclaurin series for exp(t)
# (the multiplication by 2**n just amounts to shifting the exponent).

# For very high precision use the newton method to compute exp from
# log_agm; for |x| very large or very small use
# exp(x + m) = exp(x) * e**m,  m = int(n * math.log(2))

# Input: x * 2**prec
# Output: exp(x) * 2**(prec + r)
def exp_series(x, prec, r):
    x >>= r
    # 1 + x + x^2/2! + x^3/3! + x^4/4! + ... =
    # (1 + x^2/2! + ...) + x * (1 + x^2/3! + ...)
    s0 = s1 = (MP_ONE << prec)
    k = 2
    a = x2 = (x * x) >> prec
    while a:
        a = a // k
        s0 += a
        k += 1
        a = a // k
        s1 += a
        a = (a * x2) >> prec
        k += 1
    # Calculate s**(2**r) by repeated squaring
    s1 = (s1 * x) >> prec
    s = s0 + s1
    while r:
        s = (s*s) >> prec
        r -= 1
    return s

def exp_series2(x, prec, r):
    x >>= r
    sign = 0
    if x < 0:
        sign = 1
        x = -x
    x2 = (x*x) >> prec
    if prec < 1500:
        s1 = a = x
        k = 3
        while a:
            a = ((a * x2) >> prec) // (k*(k-1))
            s1 += a
            k += 2
    else:
        # use Smith's method:
        # reduce the number of multiplication summing concurrently J series
        # J=4
        # Sinh(x) =
        #   (x + x^9/9! + ...) + x^2 * (x/3! + x^9/11! + ...) +
        #   x^4 * (x/5! + x^9/13! + ...) + x^6 * (x/7! + x^9/15! + ...)
        J = 4
        ax = [MP_ONE << prec, x2]
        px = x2
        asum = [x, x//6]
        fact = 6
        k = 4
        for j in range(2, J):
            px = (px * x2) >> prec
            ax.append(px)
            fact *= k*(k+1)
            asum.append(x//fact)
            k += 2
        lx = (ax[-1]*x2) >> prec
        p = asum[-1]
        while p:
            p = (p * lx) >> prec
            for j in range(J):
                p = p//(k*(k+1))
                asum[j] += p
                k += 2
        s1 = 0
        for i in range(1, J):
            s1 += ax[i]*asum[i]
        s1 = asum[0] + (s1 >> prec)
    c1 = isqrt_fast((s1*s1) + (MP_ONE<<(2*prec)))
    if sign:
        s = c1 - s1
    else:
        s = c1 + s1
    # Calculate s**(2**r) by repeated squaring
    while r:
        s = (s*s) >> prec
        r -= 1
    return s

# use the fourth order newton method, with step
# r = r + r * (h + h^2/2 + h^3/6 + h$/24)
# at each step the precision is quadrupled.

def exp_newton(x, prec):
    extra = 10
    r = mpf_exp(x, 60)
    start = 50
    prevp = start
    for p in giant_steps(start, prec+extra, 4):
        h = mpf_sub(x, mpf_log(r, p), p)
        h2 = mpf_mul(h, h, p)
        h3 = mpf_mul(h2, h, p)
        h4 = mpf_mul(h2, h2, p)
        t = mpf_add(h, mpf_shift(h2, -1), p)
        t = mpf_add(t, mpf_div(h3, from_int(6, p), p), p)
        t = mpf_add(t, mpf_div(h4, from_int(24, p), p), p)
        t = mpf_mul(r, t, p)
        r = mpf_add(r, t, p)
    return r

# for precision larger than this limit, for x > 1, use the newton method
LIM_EXP_SERIES2 = 10000
# when the newton method is used, if x has mag=exp+bc larger than LIM_MAG
# shift it
LIM_MAG = 5

# table of values to determine if exp_series2 or exp_newton is faster,
# determined with benchmarking on a PC, with gmpy
ns_exp = [8,9,10,11,12,13,33,66,83,99,132,166,199,232,265,298,332,664]
precs_exp = [43000, 63000, 64000, 64000, 65000, 66000, 72000, 82000, 99000,
   115000, 148000, 190000, 218000, 307000, 363000, 528000, 594000, 1650000]


def mpf_exp(x, prec, rnd=round_fast):
    sign, man, exp, bc = x
    if not man:
        if not exp:
            return fone
        if x == fninf:
            return fzero
        return x
    mag = bc+exp
    # Fast handling e**n. TODO: the best cutoff depends on both the
    # size of n and the precision.
    if prec > 600 and exp >= 0:
        e = mpf_e(prec+10+int(1.45*mag))
        return mpf_pow_int(e, (-1)**sign *(man<<exp), prec, rnd)
    if mag < -prec-10:
        return mpf_perturb(fone, sign, prec, rnd)
    # extra precision needs to be similar in magnitude to log_2(|x|)
    # for the modulo reduction, plus r for the error from squaring r times
    wp = prec + max(0, mag)
    if wp < 300:
        r = int(2*wp**0.4)
        if mag < 0:
            r = max(1, r + mag)
        wp += r + 20
        t = to_fixed(x, wp)
        # abs(x) > 1?
        if mag > 1:
            lg2 = ln2_fixed(wp)
            n, t = divmod(t, lg2)
        else:
            n = 0
        man = exp_series(t, wp, r)
    else:
        use_newton = False
        # put a bound on exp to avoid infinite recursion in exp_newton
        # TODO find a good bound
        if wp > LIM_EXP_SERIES2 and exp < 1000:
            if mag > 0:
                use_newton = True
            elif mag <= 0 and -mag <= ns_exp[-1]:
                i = bisect(ns_exp, -mag-1)
                if i < len(ns_exp):
                    wp0 = precs_exp[i]
                    if wp > wp0:
                        use_newton = True

        if not use_newton:
            r = int(0.7 * wp**0.5)
            if mag < 0:
                r = max(1, r + mag)
            wp += r + 20
            t = to_fixed(x, wp)
            if mag > 1:
                lg2 = ln2_fixed(wp)
                n, t = divmod(t, lg2)
            else:
                n = 0
            man = exp_series2(t, wp, r)
        else:
            # if x is very small or very large use
            # exp(x + m) = exp(x) * e**m
            if mag > LIM_MAG:
                wp += mag*10 + 100
                n = int(mag * math.log(2)) + 1
                x = mpf_sub(x, from_int(n, wp), wp)
            elif mag <= 0:
                wp += -mag*10 + 100
                if mag < 0:
                    n = int(-mag * math.log(2)) + 1
                    x = mpf_add(x, from_int(n, wp), wp)
            res = exp_newton(x, wp)
            sign, man, exp, bc = res
            if mag < 0:
                t = mpf_pow_int(mpf_e(wp), n, wp)
                res = mpf_div(res, t, wp)
                sign, man, exp, bc = res
            if mag > LIM_MAG:
                t = mpf_pow_int(mpf_e(wp), n, wp)
                res = mpf_mul(res, t, wp)
                sign, man, exp, bc = res
            return normalize(sign, man, exp, bc, prec, rnd)
    bc = bitcount(man)
    return normalize(0, man, int(-wp+n), bc, prec, rnd)


#----------------------------------------------------------------------------#
#                                                                            #
#                          Trigonometric functions                           #
#                                                                            #
#----------------------------------------------------------------------------#

def sin_taylor(x, prec):
    x = MP_BASE(x)
    x2 = (x*x) >> prec
    s = a = x
    k = 3
    while a:
        a = ((a * x2) >> prec) // (k*(1-k))
        s += a
        k += 2
    return s

def cos_taylor(x, prec):
    x = MP_BASE(x)
    x2 = (x*x) >> prec
    a = c = (MP_ONE<<prec)
    k = 2
    while a:
        a = ((a * x2) >> prec) // (k*(1-k))
        c += a
        k += 2
    return c

# Input: x * 2**prec
# Output: c * 2**(prec + r), s * 2**(prec + r)
def expi_series(x, prec, r):
    x >>= r
    one = MP_ONE << prec
    x2 = (x*x) >> prec
    s = x
    a = x
    k = 2
    while a:
        a = ((a * x2) >> prec) // (-k*(k+1))
        s += a
        k += 2
    c = isqrt_fast((MP_ONE<<(2*prec)) - (s*s))
    # Calculate (c + j*s)**(2**r) by repeated squaring
    for j in range(r):
        c, s =  (c*c-s*s) >> prec, (2*c*s ) >> prec
    return c, s

def reduce_angle(x, prec):
    """
    Let x be a nonzero, finite mpf value defining angle (measured in
    radians). Then reduce_trig(x, prec) returns (y, swaps, n) where:

      y = (man, wp) is the reduced angle as a scaled fixed-point
        number with precision wp, i.e. a floating-point number with
        exponent -wp. The mantissa is positive and has width ~equal
        to the input prec.

      swaps = (swap_cos_sin, cos_sign, sin_sign)
        Flags indicating the swaps that need to be applied
        to (cos(y), sin(y)) to obtain (cos(x), sin(x))

      n is an integer giving the original quadrant of x

    Calculation of the quadrant
    ===========================

    The integer n indices the quadrant of x. That is:

        ...
        -pi     <   x  < -pi/2     n = -2
        -pi/2   <   x  <  0        n = -1
        0       <   x  <  pi/2     n = 0
        pi/2    <   x  <  pi       n = 1
        pi      <   x  <  3*pi/2   n = 2
        3*pi/2  <   x  <  2*pi     n = 3
        2*pi    <   x  <  5*pi/2   n = 4
        ...

    Note that n does not wrap around. A quadrant index normalized to
    lie in [0, 1, 2, 3] can be found easily later on by computing
    n % 4. Keeping the extended information in n is crucial for
    interval arithmetic, as it allows one to distinguish between
    whether two points of a sine wave lie next to each other on
    a monotonic segment or are actually separated by a full
    period (or several periods).

    Note also that because is x is guaranteed to be rational, and
    all roots of the sine/cosine are irrational, all inequalities are
    strict. That is, we can always compute the correct quadrant.
    Care is required to do ensure that this is done right.

    Swaps
    =====

    The number y is a reduction of x to the first quadrant. This is
    essentially x mod pi/2. In fact, we reduce y further, to the first
    octant, by computing pi/2-x if x > pi/4.

    Due to the translation and mirror symmetries of trigonometric
    functions, this allows us to compute sin(x) or cos(x) by computing
    +/-sin(y) or +/-cos(y). The point, of course, is that if x
    is large, the Taylor series for y converges much more quickly
    than the one for x.

    """
    sign, man, exp, bc = x
    magnitude = exp + bc

    if not man:
        return (0, 0), (0, 0, 0), 0

    # Here we have abs(x) < 0.5. In this case no reduction is necessary.
    # TODO: could also handle abs(x) < 1
    if magnitude < 0:
        # Quadrant is 0 or -1
        n = -sign
        swaps = (0, 0, sign)
        fixed_exp = exp + bc - prec
        delta = fixed_exp - exp
        if delta < 0:
            man <<= (-delta)
        elif delta > 0:
            man >>= delta
        y = (man, -fixed_exp)
        return y, swaps, n

    i = 0
    while 1:
        cancellation_prec = 20 * 2**i
        wp = prec + abs(magnitude) + cancellation_prec
        pi1 = pi_fixed(wp)
        pi2 = pi1 >> 1
        pi4 = pi1 >> 2
        # Find nearest multiple
        n, y = divmod(to_fixed(x, wp), pi2)
        # Interchange cos/sin ?
        if y > pi4:
            swap_cos_sin = 1
            y = pi2 - y
        else:
            swap_cos_sin = 0
        # Now, the catch is that x might be extremely close to a
        # multiple of pi/2. This means accuracy is lost, and we may
        # even end up in the wrong quadrant, which is bad news
        # for interval arithmetic. This effect manifests by the
        # fixed-point value of y becoming small.  This is easy to check for.
        if y >> (prec + magnitude - 10):
            n = int(n)
            swaps = swap_table[swap_cos_sin^(n%2)][n%4]
            return (y>>magnitude, wp-magnitude), swaps, n
        i += 1

swap_table = ((0,0,0),(0,1,0),(0,1,1),(0,0,1)), ((1,0,0),(1,1,0),(1,1,1),(1,0,1))

def calc_cos_sin(which, y, swaps, prec, cos_rnd, sin_rnd):
    """
    Simultaneous computation of cos and sin (internal function).
    """
    y, wp = y
    swap_cos_sin, cos_sign, sin_sign = swaps

    if swap_cos_sin:
        which_compute = -which
    else:
        which_compute = which

    # XXX: assumes no swaps
    if not y:
        return fone, fzero

    # Tiny nonzero argument
    if wp > prec*2 + 30:
        y = from_man_exp(y, -wp)

        if swap_cos_sin:
            cos_rnd, sin_rnd = sin_rnd, cos_rnd
            cos_sign, sin_sign = sin_sign, cos_sign

        if cos_sign: cos = mpf_perturb(fnone, 0, prec, cos_rnd)
        else:        cos = mpf_perturb(fone, 1, prec, cos_rnd)
        if sin_sign: sin = mpf_perturb(mpf_neg(y), 0, prec, sin_rnd)
        else:        sin = mpf_perturb(y, 1, prec, sin_rnd)

        if swap_cos_sin:
            cos, sin = sin, cos
        return cos, sin

    # Use standard Taylor series
    if prec < 600:
        if which_compute == 0:
            sin = sin_taylor(y, wp)
            # only need to evaluate one of the series
            cos = isqrt_fast((MP_ONE<<(2*wp)) - sin*sin)
        elif which_compute == 1:
            sin = 0
            cos = cos_taylor(y, wp)
        elif which_compute == -1:
            sin = sin_taylor(y, wp)
            cos = 0
    # Use exp(i*x) with Brent's trick
    else:
        r = int(0.137 * prec**0.579)
        ep = r+20
        cos, sin = expi_series(y<<ep, wp+ep, r)
        cos >>= ep
        sin >>= ep

    if swap_cos_sin:
        cos, sin = sin, cos

    if cos_rnd is not round_nearest:
        # Round and set correct signs
        # XXX: this logic needs a second look
        ONE = MP_ONE << wp
        if cos_sign:
            cos += (-1)**(cos_rnd in (round_ceiling, round_down))
            cos = min(ONE, cos)
        else:
            cos += (-1)**(cos_rnd in (round_ceiling, round_up))
            cos = min(ONE, cos)
        if sin_sign:
            sin += (-1)**(sin_rnd in (round_ceiling, round_down))
            sin = min(ONE, sin)
        else:
            sin += (-1)**(sin_rnd in (round_ceiling, round_up))
            sin = min(ONE, sin)

    if which != -1:
        cos = normalize(cos_sign, cos, -wp, bitcount(cos), prec, cos_rnd)
    if which != 1:
        sin = normalize(sin_sign, sin, -wp, bitcount(sin), prec, sin_rnd)

    return cos, sin

def cos_sin(x, prec, rnd=round_fast, which=0):
    """
    Computes (cos(x), sin(x)). The parameter 'which' can disable
    evaluation of either cos or sin:

        0 -- return (cos(x), sin(x), n)
        1 -- return (cos(x), -,      n)
       -1 -- return (-,      sin(x), n)

    If only one function is wanted, this is slightly
    faster at low precision.
    """
    sign, man, exp, bc = x
    # Exact (or special) cases
    if not man:
        if exp:
            return (fnan, fnan)
        else:
            return (fone, fzero)
    y, swaps, n = reduce_angle(x, prec+10)
    return calc_cos_sin(which, y, swaps, prec, rnd, rnd)

def mpf_cos(x, prec, rnd=round_fast):
    return cos_sin(x, prec, rnd, 1)[0]

def mpf_sin(x, prec, rnd=round_fast):
    return cos_sin(x, prec, rnd, -1)[1]

def mpf_tan(x, prec, rnd=round_fast):
    c, s = cos_sin(x, prec+20)
    return mpf_div(s, c, prec, rnd)

# Accurate computation of cos(pi*x) and sin(pi*x) is needed by
# reflection formulas for gamma, polygamma, zeta, etc

def mpf_cos_sin_pi(x, prec, rnd=round_fast):
    """Accurate computation of (cos(pi*x), sin(pi*x))
    for x close to an integer"""
    sign, man, exp, bc = x
    if not man:
        return cos_sin(x, prec, rnd)
    # Exactly an integer or half-integer?
    if exp >= -1:
        if exp == -1:
            c = fzero
            s = (fone, fnone)[bool(man & 2) ^ sign]
        elif exp == 0:
            c, s = (fnone, fzero)
        else:
            c, s = (fone, fzero)
        return c, s
    # Close to 0 ?
    size = exp + bc
    if size < -(prec+5):
        c = mpf_perturb(fone, 1, prec, rnd)
        s = mpf_perturb(mpf_mul(x, mpf_pi(prec)), sign, prec, rnd)
        return c, s
    if sign:
        man = -man
    # Subtract nearest half-integer (= modulo pi/2)
    nhint = ((man >> (-exp-2)) + 1) >> 1
    man = man - (nhint << (-exp-1))
    x = from_man_exp(man, exp, prec)
    x = mpf_mul(x, mpf_pi(prec), prec)
    # XXX: with some more work, could call calc_cos_sin,
    # to save some time and to get rounding right
    case = nhint % 4
    if case == 0:
        c, s = cos_sin(x, prec, rnd)
    elif case == 1:
        s, c = cos_sin(x, prec, rnd)
        c = mpf_neg(c)
    elif case == 2:
        c, s = cos_sin(x, prec, rnd)
        c = mpf_neg(c)
        s = mpf_neg(s)
    else:
        s, c = cos_sin(x, prec, rnd)
        s = mpf_neg(s)
    return c, s

def mpf_cos_pi(x, prec, rnd=round_fast):
    return mpf_cos_sin_pi(x, prec, rnd)[0]

def mpf_sin_pi(x, prec, rnd=round_fast):
    return mpf_cos_sin_pi(x, prec, rnd)[1]


#----------------------------------------------------------------------
# Hyperbolic functions
#

def sinh_taylor(x, prec):
    x = MP_BASE(x)
    x2 = (x*x) >> prec
    s = a = x
    k = 3
    while a:
        a = ((a * x2) >> prec) // (k*(k-1))
        s += a
        k += 2
    return s

def cosh_sinh(x, prec, rnd=round_fast, tanh=0):
    """Simultaneously compute (cosh(x), sinh(x)) for real x"""
    sign, man, exp, bc = x
    if (not man) and exp:
        if tanh:
            if x == finf: return fone
            if x == fninf: return fnone
            return fnan
        if x == finf: return (finf, finf)
        if x == fninf: return (finf, fninf)
        return fnan, fnan

    if sign:
        man = -man

    mag = exp + bc
    prec2 = prec + 20

    if mag < -3:
        # Extremely close to 0, sinh(x) ~= x and cosh(x) ~= 1
        if mag < -prec-2:
            if tanh:
                return mpf_perturb(x, 1-sign, prec, rnd)
            cosh = mpf_perturb(fone, 0, prec, rnd)
            sinh = mpf_perturb(x, sign, prec, rnd)
            return cosh, sinh

        # Avoid cancellation when computing sinh
        # TODO: might be faster to use sinh series directly
        prec2 += (-mag) + 4

    # In the general case, we use
    #    cosh(x) = (exp(x) + exp(-x))/2
    #    sinh(x) = (exp(x) - exp(-x))/2
    # and note that the exponential only needs to be computed once.
    ep = mpf_exp(x, prec2)
    em = mpf_div(fone, ep, prec2)
    if tanh:
        ch = mpf_add(ep, em, prec2, rnd)
        sh = mpf_sub(ep, em, prec2, rnd)
        return mpf_div(sh, ch, prec, rnd)
    else:
        ch = mpf_shift(mpf_add(ep, em, prec, rnd), -1)
        sh = mpf_shift(mpf_sub(ep, em, prec, rnd), -1)
        return ch, sh

def mpf_cosh(x, prec, rnd=round_fast):
    """Compute cosh(x) for a real argument x"""
    return cosh_sinh(x, prec, rnd)[0]

def mpf_sinh(x, prec, rnd=round_fast):
    """Compute sinh(x) for a real argument x"""
    return cosh_sinh(x, prec, rnd)[1]

def mpf_tanh(x, prec, rnd=round_fast):
    """Compute tanh(x) for a real argument x"""
    return cosh_sinh(x, prec, rnd, tanh=1)


#----------------------------------------------------------------------
# Inverse tangent
#

def atan_newton(x, prec):
    if prec >= 100:
        r = math.atan((x>>(prec-53))/2.0**53)
    else:
        r = math.atan(x/2.0**prec)
    prevp = 50
    r = int(r * 2.0**53) >> (53-prevp)
    extra_p = 100
    for p in giant_steps(prevp, prec):
        s = int(0.137 * p**0.579)
        p += s + 50
        r = r << (p-prevp)
        cos, sin = expi_series(r, p, s)
        tan = (sin << p) // cos
        a = ((tan - rshift(x, prec-p)) << p) // ((MP_ONE<<p) + ((tan**2)>>p))
        r = r - a
        prevp = p
    return rshift(r, prevp-prec)


ATAN_TAYLOR_PREC = 3000
ATAN_TAYLOR_SHIFT = 7   # steps of size 2^-N

atan_taylor_cache = {}

def atan_taylor_get_cached(n, prec):
    # Taylor series with caching wins up to huge precisions
    # To avoid unnecessary precomputation at low precision, we
    # do it in steps
    # Round to next power of 2
    prec2 = (1<<(bitcount(prec-1))) + 20
    dprec = prec2 - prec
    if (n, prec2) in atan_taylor_cache:
        a, atan_a = atan_taylor_cache[n, prec2]
    else:
        a = n << (prec2 - ATAN_TAYLOR_SHIFT)
        atan_a = atan_newton(a, prec2)
        atan_taylor_cache[n, prec2] = (a, atan_a)
    return (a >> dprec), (atan_a >> dprec)

def atan_taylor(x, prec):
    n = (x >> (prec-ATAN_TAYLOR_SHIFT))
    a, atan_a = atan_taylor_get_cached(n, prec)
    d = x - a
    s0 = v = (d << prec) // ((a**2 >> prec) + (a*d >> prec) + (MP_ONE << prec))
    v2 = (v**2 >> prec)
    v4 = (v2 * v2) >> prec
    s1 = v//3
    v = (v * v4) >> prec
    k = 5
    while v:
        s0 += v // k
        k += 2
        s1 += v // k
        v = (v * v4) >> prec
        k += 2
    s1 = (s1 * v2) >> prec
    s = s0 - s1
    return atan_a + s

def atan_inf(sign, prec, rnd):
    if not sign:
        return mpf_shift(mpf_pi(prec, rnd), -1)
    return mpf_neg(mpf_shift(mpf_pi(prec, negative_rnd[rnd]), -1))

def mpf_atan(x, prec, rnd=round_fast):
    sign, man, exp, bc = x
    if not man:
        if x == fzero: return fzero
        if x == finf: return atan_inf(0, prec, rnd)
        if x == fninf: return atan_inf(1, prec, rnd)
        return fnan
    mag = exp + bc
    # Essentially infinity
    if mag > prec+20:
        return atan_inf(sign, prec, rnd)
    # Essentially ~ x
    if -mag > prec+20:
        return mpf_perturb(x, 1-sign, prec, rnd)
    wp = prec + 30 + abs(mag)
    # For large x, use atan(x) = pi/2 - atan(1/x)
    if mag >= 2:
        x = mpf_rdiv_int(1, x, wp)
        reciprocal = True
    else:
        reciprocal = False
    t = to_fixed(x, wp)
    if sign:
        t = -t
    if wp < ATAN_TAYLOR_PREC:
        a = atan_taylor(t, wp)
    else:
        a = atan_newton(t, wp)
    if reciprocal:
        a = ((pi_fixed(wp)>>1)+1) - a
    if sign:
        a = -a
    return from_man_exp(a, -wp, prec, rnd)

# TODO: cleanup the special cases
def mpf_atan2(y, x, prec, rnd=round_fast):
    xsign, xman, xexp, xbc = x
    ysign, yman, yexp, ybc = y
    if not yman:
        if y == fzero and x != fnan:
            if mpf_sign(x) >= 0:
                return fzero
            return mpf_pi(prec, rnd)
        if y in (finf, fninf):
            if x in (finf, fninf):
                return fnan
            # pi/2
            if y == finf:
                return mpf_shift(mpf_pi(prec, rnd), -1)
            # -pi/2
            return mpf_neg(mpf_shift(mpf_pi(prec, negative_rnd[rnd]), -1))
        return fnan
    if ysign:
        return mpf_neg(mpf_atan2(mpf_neg(y), x, prec, negative_rnd[rnd]))
    if not xman:
        if x == fnan:
            return fnan
        if x == finf:
            return fzero
        if x == fninf:
            return mpf_pi(prec, rnd)
        if y == fzero:
            return fzero
        return mpf_shift(mpf_pi(prec, rnd), -1)
    tquo = mpf_atan(mpf_div(y, x, prec+4), prec+4)
    if xsign:
        return mpf_add(mpf_pi(prec+4), tquo, prec, rnd)
    else:
        return mpf_pos(tquo, prec, rnd)

def mpf_asin(x, prec, rnd=round_fast):
  sign, man, exp, bc = x
  if bc+exp > 0 and x not in (fone, fnone):
      raise ComplexResult("asin(x) is real only for -1 <= x <= 1")
  # asin(x) = 2*atan(x/(1+sqrt(1-x**2)))
  wp = prec + 15
  a = mpf_mul(x, x)
  b = mpf_add(fone, mpf_sqrt(mpf_sub(fone, a, wp), wp), wp)
  c = mpf_div(x, b, wp)
  return mpf_shift(mpf_atan(c, prec, rnd), 1)

def mpf_acos(x, prec, rnd=round_fast):
    # acos(x) = 2*atan(sqrt(1-x**2)/(1+x))
    sign, man, exp, bc = x
    if bc + exp > 0:
        if x not in (fone, fnone):
            raise ComplexResult("acos(x) is real only for -1 <= x <= 1")
        if x == fnone:
            return mpf_pi(prec, rnd)
    wp = prec + 15
    a = mpf_mul(x, x)
    b = mpf_sqrt(mpf_sub(fone, a, wp), wp)
    c = mpf_div(b, mpf_add(fone, x, wp), wp)
    return mpf_shift(mpf_atan(c, prec, rnd), 1)

def mpf_asinh(x, prec, rnd=round_fast):
    wp = prec + 20
    sign, man, exp, bc = x
    mag = exp+bc
    if mag < -8:
        if mag < -wp:
            return mpf_perturb(x, 1-sign, prec, rnd)
        wp += (-mag)
    # asinh(x) = log(x+sqrt(x**2+1))
    # use reflection symmetry to avoid cancellation
    q = mpf_sqrt(mpf_add(mpf_mul(x, x), fone, wp), wp)
    q = mpf_add(mpf_abs(x), q, wp)
    if sign:
        return mpf_neg(mpf_log(q, prec, negative_rnd[rnd]))
    else:
        return mpf_log(q, prec, rnd)

def mpf_acosh(x, prec, rnd=round_fast):
    # acosh(x) = log(x+sqrt(x**2-1))
    wp = prec + 15
    if mpf_cmp(x, fone) == -1:
        raise ComplexResult("acosh(x) is real only for x >= 1")
    q = mpf_sqrt(mpf_add(mpf_mul(x,x), fnone, wp), wp)
    return mpf_log(mpf_add(x, q, wp), prec, rnd)

def mpf_atanh(x, prec, rnd=round_fast):
    # atanh(x) = log((1+x)/(1-x))/2
    sign, man, exp, bc = x
    if (not man) and exp:
        if x in (fzero, fnan):
            return x
        raise ComplexResult("atanh(x) is real only for -1 <= x <= 1")
    mag = bc + exp
    if mag > 0:
        if mag == 1 and man == 1:
            return [finf, fninf][sign]
        raise ComplexResult("atanh(x) is real only for -1 <= x <= 1")
    wp = prec + 15
    if mag < -8:
        if mag < -wp:
            return mpf_perturb(x, sign, prec, rnd)
        wp += (-mag)
    a = mpf_add(x, fone, wp)
    b = mpf_sub(fone, x, wp)
    return mpf_shift(mpf_log(mpf_div(a, b, wp), prec, rnd), -1)

def mpf_fibonacci(x, prec, rnd=round_fast):
    sign, man, exp, bc = x
    if not man:
        if x == fninf:
            return fnan
        return x
    # F(2^n) ~= 2^(2^n)
    size = abs(exp+bc)
    if exp >= 0:
        # Exact
        if size < 10 or size <= bitcount(prec):
            return from_int(ifib(to_int(x)), prec, rnd)
    # Use the modified Binet formula
    wp = prec + size + 20
    a = mpf_phi(wp)
    b = mpf_add(mpf_shift(a, 1), fnone, wp)
    u = mpf_pow(a, x, wp)
    v = mpf_cos_pi(x, wp)
    v = mpf_div(v, u, wp)
    u = mpf_sub(u, v, wp)
    u = mpf_div(u, b, prec, rnd)
    return u
