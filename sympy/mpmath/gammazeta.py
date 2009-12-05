"""
-----------------------------------------------------------------------
This module implements gamma- and zeta-related functions:

* Bernoulli numbers
* Factorials
* The gamma function
* Polygamma functions
* Harmonic numbers
* The Riemann zeta function
* Constants related to these functions

-----------------------------------------------------------------------
"""

import math

from settings import (\
    MP_BASE, MP_ZERO, MP_ONE, MP_THREE, MODE, gmpy,
    round_floor, round_ceiling, round_down, round_up,
    round_nearest, round_fast
)

from libintmath import list_primes, int_fac, moebius

from libmpf import (\
    lshift, sqrt_fixed,
    fzero, fone, fnone, fhalf, ftwo, finf, fninf, fnan,
    from_int, to_int, to_fixed, from_man_exp, from_rational,
    mpf_pos, mpf_neg, mpf_abs, mpf_add, mpf_sub,
    mpf_mul, mpf_mul_int, mpf_div, mpf_sqrt, mpf_pow_int,
    mpf_rdiv_int,
    mpf_perturb, mpf_le, mpf_lt, mpf_shift,
    negative_rnd, reciprocal_rnd,
)

from libelefun import (\
    constant_memo,
    def_mpf_constant,
    mpf_pi, pi_fixed, ln2_fixed, log_int_fixed, mpf_ln2,
    mpf_exp, mpf_log, mpf_pow, mpf_cosh,
    cos_sin, cosh_sinh, mpf_cos_sin_pi, mpf_cos_pi, mpf_sin_pi,
)

from libmpc import (\
    mpc_zero, mpc_one, mpc_half, mpc_two,
    mpc_abs, mpc_shift, mpc_pos,
    mpc_add, mpc_sub, mpc_mul, mpc_div,
    mpc_add_mpf, mpc_mul_mpf, mpc_div_mpf, mpc_mpf_div,
    mpc_mul_int, mpc_pow_int,
    mpc_log, mpc_exp, mpc_pow,
    mpc_cos_pi, mpc_sin_pi,
    mpc_reciprocal, mpc_square
)

# Catalan's constant is computed using Lupas's rapidly convergent series
# (listed on http://mathworld.wolfram.com/CatalansConstant.html)
#            oo
#            ___       n-1  8n     2                   3    2
#        1  \      (-1)    2   (40n  - 24n + 3) [(2n)!] (n!)
#  K =  ---  )     -----------------------------------------
#       64  /___               3               2
#                             n  (2n-1) [(4n)!]
#           n = 1

@constant_memo
def catalan_fixed(prec):
    prec = prec + 20
    a = one = MP_ONE << prec
    s, t, n = 0, 1, 1
    while t:
        a *= 32 * n**3 * (2*n-1)
        a //= (3-16*n+16*n**2)**2
        t = a * (-1)**(n-1) * (40*n**2-24*n+3) // (n**3 * (2*n-1))
        s += t
        n += 1
    return s >> (20 + 6)

# Khinchin's constant is relatively difficult to compute. Here
# we use the rational zeta series

#                    oo                2*n-1
#                   ___                ___
#                   \   ` zeta(2*n)-1  \   ` (-1)^(k+1)
#  log(K)*log(2) =   )    ------------  )    ----------
#                   /___.      n       /___.      k
#                   n = 1              k = 1

# which adds half a digit per term. The essential trick for achieving
# reasonable efficiency is to recycle both the values of the zeta
# function (essentially Bernoulli numbers) and the partial terms of
# the inner sum.

# An alternative might be to use K = 2*exp[1/log(2) X] where

#      / 1     1       [ pi*x*(1-x^2) ]
#  X = |    ------ log [ ------------ ].
#      / 0  x(1+x)     [  sin(pi*x)   ]

# and integrate numerically. In practice, this seems to be slightly
# slower than the zeta series at high precision.

@constant_memo
def khinchin_fixed(prec):
    wp = int(prec + prec**0.5 + 15)
    s = MP_ZERO
    fac = from_int(4)
    t = ONE = MP_ONE << wp
    pi = mpf_pi(wp)
    pipow = twopi2 = mpf_shift(mpf_mul(pi, pi, wp), 2)
    n = 1
    while 1:
        zeta2n = mpf_abs(mpf_bernoulli(2*n, wp))
        zeta2n = mpf_mul(zeta2n, pipow, wp)
        zeta2n = mpf_div(zeta2n, fac, wp)
        zeta2n = to_fixed(zeta2n, wp)
        term = (((zeta2n - ONE) * t) // n) >> wp
        if term < 100:
            break
        #if not n % 10:
        #    print n, math.log(int(abs(term)))
        s += term
        t += ONE//(2*n+1) - ONE//(2*n)
        n += 1
        fac = mpf_mul_int(fac, (2*n)*(2*n-1), wp)
        pipow = mpf_mul(pipow, twopi2, wp)
    s = (s << wp) // ln2_fixed(wp)
    K = mpf_exp(from_man_exp(s, -wp), wp)
    K = to_fixed(K, prec)
    return K


# Glaisher's constant is defined as A = exp(1/2 - zeta'(-1)).
# One way to compute it would be to perform direct numerical
# differentiation, but computing arbitrary Riemann zeta function
# values at high precision is expensive. We instead use the formula

#     A = exp((6 (-zeta'(2))/pi^2 + log 2 pi + gamma)/12)

# and compute zeta'(2) from the series representation

#              oo
#              ___
#             \     log k
#  -zeta'(2) = )    -----
#             /___     2
#                    k
#            k = 2

# This series converges exceptionally slowly, but can be accelerated
# using Euler-Maclaurin formula. The important insight is that the
# E-M integral can be done in closed form and that the high order
# are given by

#    n  /       \
#   d   | log x |   a + b log x
#   --- | ----- | = -----------
#     n |   2   |      2 + n
#   dx  \  x    /     x

# where a and b are integers given by a simple recurrence. Note
# that just one logarithm is needed. However, lots of integer
# logarithms are required for the initial summation.

# This algorithm could possibly be turned into a faster algorithm
# for general evaluation of zeta(s) or zeta'(s); this should be
# looked into.

@constant_memo
def glaisher_fixed(prec):
    wp = prec + 30
    # Number of direct terms to sum before applying the Euler-Maclaurin
    # formula to the tail. TODO: choose more intelligently
    N = int(0.33*prec + 5)
    ONE = MP_ONE << wp
    # Euler-Maclaurin, step 1: sum log(k)/k**2 for k from 2 to N-1
    s = MP_ZERO
    for k in range(2, N):
        #print k, N
        s += log_int_fixed(k, wp) // k**2
    logN = log_int_fixed(N, wp)
    #logN = to_fixed(mpf_log(from_int(N), wp+20), wp)
    # E-M step 2: integral of log(x)/x**2 from N to inf
    s += (ONE + logN) // N
    # E-M step 3: endpoint correction term f(N)/2
    s += logN // (N**2 * 2)
    # E-M step 4: the series of derivatives
    pN = N**3
    a = 1
    b = -2
    j = 3
    fac = from_int(2)
    k = 1
    while 1:
        # D(2*k-1) * B(2*k) / fac(2*k) [D(n) = nth derivative]
        D = ((a << wp) + b*logN) // pN
        D = from_man_exp(D, -wp)
        B = mpf_bernoulli(2*k, wp)
        term = mpf_mul(B, D, wp)
        term = mpf_div(term, fac, wp)
        term = to_fixed(term, wp)
        if abs(term) < 100:
            break
        #if not k % 10:
        #    print k, math.log(int(abs(term)), 10)
        s -= term
        # Advance derivative twice
        a, b, pN, j = b-a*j, -j*b, pN*N, j+1
        a, b, pN, j = b-a*j, -j*b, pN*N, j+1
        k += 1
        fac = mpf_mul_int(fac, (2*k)*(2*k-1), wp)
    # A = exp((6*s/pi**2 + log(2*pi) + euler)/12)
    pi = pi_fixed(wp)
    s *= 6
    s = (s << wp) // (pi**2 >> wp)
    s += euler_fixed(wp)
    s += to_fixed(mpf_log(from_man_exp(2*pi, -wp), wp), wp)
    s //= 12
    A = mpf_exp(from_man_exp(s, -wp), wp)
    return to_fixed(A, prec)

# Apery's constant can be computed using the very rapidly convergent
# series
#              oo
#              ___              2                      10
#             \         n  205 n  + 250 n + 77     (n!)
#  zeta(3) =   )    (-1)   -------------------  ----------
#             /___               64                      5
#             n = 0                             ((2n+1)!)

@constant_memo
def apery_fixed(prec):
    prec += 20
    d = MP_ONE << prec
    term = MP_BASE(77) << prec
    n = 1
    s = MP_ZERO
    while term:
        s += term
        d *= (n**10)
        d //= (((2*n+1)**5) * (2*n)**5)
        term = (-1)**n * (205*(n**2) + 250*n + 77) * d
        n += 1
    return s >> (20 + 6)

"""
Euler's constant (gamma) is computed using the Brent-McMillan formula,
gamma ~= I(n)/J(n) - log(n), where

   I(n) = sum_{k=0,1,2,...} (n**k / k!)**2 * H(k)
   J(n) = sum_{k=0,1,2,...} (n**k / k!)**2
   H(k) = 1 + 1/2 + 1/3 + ... + 1/k

The error is bounded by O(exp(-4n)). Choosing n to be a power
of two, 2**p, the logarithm becomes particularly easy to calculate.[1]

We use the formulation of Algorithm 3.9 in [2] to make the summation
more efficient.

Reference:
[1] Xavier Gourdon & Pascal Sebah, The Euler constant: gamma
http://numbers.computation.free.fr/Constants/Gamma/gamma.pdf

[2] Jonathan Borwein & David Bailey, Mathematics by Experiment,
A K Peters, 2003
"""

@constant_memo
def euler_fixed(prec):
    extra = 30
    prec += extra
    # choose p such that exp(-4*(2**p)) < 2**-n
    p = int(math.log((prec/4) * math.log(2), 2)) + 1
    n = 2**p
    A = U = -p*ln2_fixed(prec)
    B = V = MP_ONE << prec
    k = 1
    while 1:
        B = B*n**2//k**2
        A = (A*n**2//k + B)//k
        U += A
        V += B
        if max(abs(A), abs(B)) < 100:
            break
        k += 1
    return (U<<(prec-extra))//V

# Use zeta accelerated formulas for the Mertens and twin
# prime constants; see
# http://mathworld.wolfram.com/MertensConstant.html
# http://mathworld.wolfram.com/TwinPrimesConstant.html

@constant_memo
def mertens_fixed(prec):
    wp = prec + 20
    m = 2
    s = mpf_euler(wp)
    while 1:
        t = mpf_zeta_int(m, wp)
        if t == fone:
            break
        t = mpf_log(t, wp)
        t = mpf_mul_int(t, moebius(m), wp)
        t = mpf_div(t, from_int(m), wp)
        s = mpf_add(s, t)
        m += 1
    return to_fixed(s, prec)

@constant_memo
def twinprime_fixed(prec):
    def I(n):
        return sum(moebius(d)<<(n//d) for d in xrange(1,n+1) if not n%d)//n
    wp = 2*prec + 30
    res = fone
    primes = [from_rational(1,p,wp) for p in [2,3,5,7]]
    ppowers = [mpf_mul(p,p,wp) for p in primes]
    n = 2
    while 1:
        a = mpf_zeta_int(n, wp)
        for i in range(4):
            a = mpf_mul(a, mpf_sub(fone, ppowers[i]), wp)
            ppowers[i] = mpf_mul(ppowers[i], primes[i], wp)
        a = mpf_pow_int(a, -I(n), wp)
        if mpf_pos(a, prec+10, 'n') == fone:
            break
        #from libmpf import to_str
        #print n, to_str(mpf_sub(fone, a), 6)
        res = mpf_mul(res, a, wp)
        n += 1
    res = mpf_mul(res, from_int(3*15*35), wp)
    res = mpf_div(res, from_int(4*16*36), wp)
    return to_fixed(res, prec)


mpf_euler = def_mpf_constant(euler_fixed)
mpf_apery = def_mpf_constant(apery_fixed)
mpf_khinchin = def_mpf_constant(khinchin_fixed)
mpf_glaisher = def_mpf_constant(glaisher_fixed)
mpf_catalan = def_mpf_constant(catalan_fixed)
mpf_mertens = def_mpf_constant(mertens_fixed)
mpf_twinprime = def_mpf_constant(twinprime_fixed)


#-----------------------------------------------------------------------#
#                                                                       #
#                          Bernoulli numbers                            #
#                                                                       #
#-----------------------------------------------------------------------#

MAX_BERNOULLI_CACHE = 3000


"""
Small Bernoulli numbers and factorials are used in numerous summations,
so it is critical for speed that sequential computation is fast and that
values are cached up to a fairly high threshhold.

On the other hand, we also want to support fast computation of isolated
large numbers. Currently, no such acceleration is provided for integer
factorials (though it is for large floating-point factorials, which are
computed via gamma if the precision is low enough).

For sequential computation of Bernoulli numbers, we use Ramanujan's formula

                           / n + 3 \
  B   =  (A(n) - S(n))  /  |       |
   n                       \   n   /

where A(n) = (n+3)/3 when n = 0 or 2 (mod 6), A(n) = -(n+3)/6
when n = 4 (mod 6), and

         [n/6]
          ___
         \      /  n + 3  \
  S(n) =  )     |         | * B
         /___   \ n - 6*k /    n-6*k
         k = 1

For isolated large Bernoulli numbers, we use the Riemann zeta function
to calculate a numerical value for B_n. The von Staudt-Clausen theorem
can then be used to optionally find the exact value of the
numerator and denominator.
"""

bernoulli_cache = {}
f3 = from_int(3)
f6 = from_int(6)

def bernoulli_size(n):
    """Accurately estimate the size of B_n (even n > 2 only)"""
    lgn = math.log(n,2)
    return int(2.326 + 0.5*lgn + n*(lgn - 4.094))

BERNOULLI_PREC_CUTOFF = bernoulli_size(MAX_BERNOULLI_CACHE)

def mpf_bernoulli(n, prec, rnd=None):
    """Computation of Bernoulli numbers (numerically)"""
    if n < 2:
        if n < 0:
            raise ValueError("Bernoulli numbers only defined for n >= 0")
        if n == 0:
            return fone
        if n == 1:
            return mpf_neg(fhalf)
    # For odd n > 1, the Bernoulli numbers are zero
    if n & 1:
        return fzero
    # If precision is extremely high, we can save time by computing
    # the Bernoulli number at a lower precision that is sufficient to
    # obtain the exact fraction, round to the exact fraction, and
    # convert the fraction back to an mpf value at the original precision
    if prec > BERNOULLI_PREC_CUTOFF and prec > bernoulli_size(n)*1.1 + 1000:
        p, q = bernfrac(n)
        return from_rational(p, q, prec, rnd or round_floor)
    if n > MAX_BERNOULLI_CACHE:
        return mpf_bernoulli_huge(n, prec, rnd)
    wp = prec + 30
    # Reuse nearby precisions
    wp += 32 - (prec & 31)
    cached = bernoulli_cache.get(wp)
    if cached:
        numbers, state = cached
        if n in numbers:
            if not rnd:
                return numbers[n]
            return mpf_pos(numbers[n], prec, rnd)
        m, bin, bin1 = state
        if n - m > 10:
            return mpf_bernoulli_huge(n, prec, rnd)
    else:
        if n > 10:
            return mpf_bernoulli_huge(n, prec, rnd)
        numbers = {0:fone}
        m, bin, bin1 = state = [2, MP_BASE(10), MP_ONE]
        bernoulli_cache[wp] = (numbers, state)
    while m <= n:
        #print m
        case = m % 6
        # Accurately estimate size of B_m so we can use
        # fixed point math without using too much precision
        szbm = bernoulli_size(m)
        s = 0
        sexp = max(0, szbm)  - wp
        if m < 6:
            a = MP_ZERO
        else:
            a = bin1
        for j in xrange(1, m//6+1):
            usign, uman, uexp, ubc = u = numbers[m-6*j]
            if usign:
                uman = -uman
            s += lshift(a*uman, uexp-sexp)
            # Update inner binomial coefficient
            j6 = 6*j
            a *= ((m-5-j6)*(m-4-j6)*(m-3-j6)*(m-2-j6)*(m-1-j6)*(m-j6))
            a //= ((4+j6)*(5+j6)*(6+j6)*(7+j6)*(8+j6)*(9+j6))
        if case == 0: b = mpf_rdiv_int(m+3, f3, wp)
        if case == 2: b = mpf_rdiv_int(m+3, f3, wp)
        if case == 4: b = mpf_rdiv_int(-m-3, f6, wp)
        s = from_man_exp(s, sexp, wp)
        b = mpf_div(mpf_sub(b, s, wp), from_int(bin), wp)
        numbers[m] = b
        m += 2
        # Update outer binomial coefficient
        bin = bin * ((m+2)*(m+3)) // (m*(m-1))
        if m > 6:
            bin1 = bin1 * ((2+m)*(3+m)) // ((m-7)*(m-6))
        state[:] = [m, bin, bin1]
    return numbers[n]

def mpf_bernoulli_huge(n, prec, rnd=None):
    wp = prec + 10
    piprec = wp + int(math.log(n,2))
    v = mpf_gamma_int(n+1, wp)
    v = mpf_mul(v, mpf_zeta_int(n, wp), wp)
    v = mpf_mul(v, mpf_pow_int(mpf_pi(piprec), -n, wp))
    v = mpf_shift(v, 1-n)
    if not n & 3:
        v = mpf_neg(v)
    return mpf_pos(v, prec, rnd or round_fast)

def bernfrac(n):
    r"""
    Returns a tuple of integers `(p, q)` such that `p/q = B_n` exactly,
    where `B_n` denotes the `n`-th Bernoulli number. The fraction is
    always reduced to lowest terms. Note that for `n > 1` and `n` odd,
    `B_n = 0`, and `(0, 1)` is returned.

    **Examples**

    The first few Bernoulli numbers are exactly::

        >>> from mpmath import *
        >>> for n in range(15):
        ...     p, q = bernfrac(n)
        ...     print n, "%s/%s" % (p, q)
        ...
        0 1/1
        1 -1/2
        2 1/6
        3 0/1
        4 -1/30
        5 0/1
        6 1/42
        7 0/1
        8 -1/30
        9 0/1
        10 5/66
        11 0/1
        12 -691/2730
        13 0/1
        14 7/6

    This function works for arbitrarily large `n`::

        >>> p, q = bernfrac(10**4)
        >>> print q
        2338224387510
        >>> print len(str(p))
        27692
        >>> mp.dps = 15
        >>> print mpf(p) / q
        -9.04942396360948e+27677
        >>> print bernoulli(10**4)
        -9.04942396360948e+27677

    Note: :func:`bernoulli` computes a floating-point approximation
    directly, without computing the exact fraction first.
    This is much faster for large `n`.

    **Algorithm**

    :func:`bernfrac` works by computing the value of `B_n` numerically
    and then using the von Staudt-Clausen theorem [1] to reconstruct
    the exact fraction. For large `n`, this is significantly faster than
    computing `B_1, B_2, \ldots, B_2` recursively with exact arithmetic.
    The implementation has been tested for `n = 10^m` up to `m = 6`.

    In practice, :func:`bernfrac` appears to be about three times
    slower than the specialized program calcbn.exe [2]

    **References**

    1. MathWorld, von Staudt-Clausen Theorem:
       http://mathworld.wolfram.com/vonStaudt-ClausenTheorem.html

    2. The Bernoulli Number Page:
       http://www.bernoulli.org/

    """
    n = int(n)
    if n < 3:
        return [(1, 1), (-1, 2), (1, 6)][n]
    if n & 1:
        return (0, 1)
    q = 1
    for k in list_primes(n+1):
        if not (n % (k-1)):
            q *= k
    prec = bernoulli_size(n) + int(math.log(q,2)) + 20
    b = mpf_bernoulli(n, prec)
    p = mpf_mul(b, from_int(q))
    pint = to_int(p, round_nearest)
    return (pint, q)


#-----------------------------------------------------------------------#
#                                                                       #
#                         The gamma function                            #
#                                                                       #
#-----------------------------------------------------------------------#


"""
We compute the real factorial / gamma function using Spouge's approximation

    x! = (x+a)**(x+1/2) * exp(-x-a) * [c_0 + S(x) + eps]

where S(x) is the sum of c_k/(x+k) from k = 1 to a-1 and the coefficients
are given by

  c_0 = sqrt(2*pi)

         (-1)**(k-1)
  c_k =  ----------- (a-k)**(k-1/2) exp(-k+a),  k = 1,2,...,a-1
          (k - 1)!

As proved by Spouge, if we choose a = log(2)/log(2*pi)*n = 0.38*n, the
relative error eps is less than 2^(-n) for any x in the right complex
half-plane (assuming a > 2). In practice, it seems that a can be chosen
quite a bit lower still (30-50%); this possibility should be investigated.

For negative x, we use the reflection formula.

References:
-----------

John L. Spouge, "Computation of the gamma, digamma, and trigamma
functions", SIAM Journal on Numerical Analysis 31 (1994), no. 3, 931-944.
"""

spouge_cache = {}

def calc_spouge_coefficients(a, prec):
    wp = prec + int(a*1.4)
    c = [0] * a
    # b = exp(a-1)
    b = mpf_exp(from_int(a-1), wp)
    # e = exp(1)
    e = mpf_exp(fone, wp)
    # sqrt(2*pi)
    sq2pi = mpf_sqrt(mpf_shift(mpf_pi(wp), 1), wp)
    c[0] = to_fixed(sq2pi, prec)
    for k in xrange(1, a):
        # c[k] = ((-1)**(k-1) * (a-k)**k) * b / sqrt(a-k)
        term = mpf_mul_int(b, ((-1)**(k-1) * (a-k)**k), wp)
        term = mpf_div(term, mpf_sqrt(from_int(a-k), wp), wp)
        c[k] = to_fixed(term, prec)
        # b = b / (e * k)
        b = mpf_div(b, mpf_mul(e, from_int(k), wp), wp)
    return c

# Cached lookup of coefficients
def get_spouge_coefficients(prec):
    # This exact precision has been used before
    if prec in spouge_cache:
        return spouge_cache[prec]
    for p in spouge_cache:
        if 0.8 <= prec/float(p) < 1:
            return spouge_cache[p]
    # Here we estimate the value of a based on Spouge's inequality for
    # the relative error
    a = max(3, int(0.38*prec))  # 0.38 = log(2)/log(2*pi), ~= 1.26*n
    coefs = calc_spouge_coefficients(a, prec)
    spouge_cache[prec] = (prec, a, coefs)
    return spouge_cache[prec]

def spouge_sum_real(x, prec, a, c):
    x = to_fixed(x, prec)
    s = c[0]
    for k in xrange(1, a):
        s += (c[k] << prec) // (x + (k << prec))
    return from_man_exp(s, -prec, prec, round_floor)

# Unused: for fast computation of gamma(p/q)
def spouge_sum_rational(p, q, prec, a, c):
    s = c[0]
    for k in xrange(1, a):
        s += c[k] * q // (p+q*k)
    return from_man_exp(s, -prec, prec, round_floor)

# For a complex number a + b*I, we have
#
#        c_k          (a+k)*c_k     b * c_k
#  -------------  =   ---------  -  ------- * I
#  (a + b*I) + k          M            M
#
#                 2    2      2   2              2
#  where M = (a+k)  + b   = (a + b ) + (2*a*k + k )

def spouge_sum_complex(re, im, prec, a, c):
    re = to_fixed(re, prec)
    im = to_fixed(im, prec)
    sre, sim = c[0], 0
    mag = ((re**2)>>prec) + ((im**2)>>prec)
    for k in xrange(1, a):
        M = mag + re*(2*k) + ((k**2) << prec)
        sre += (c[k] * (re + (k << prec))) // M
        sim -= (c[k] * im) // M
    re = from_man_exp(sre, -prec, prec, round_floor)
    im = from_man_exp(sim, -prec, prec, round_floor)
    return re, im

# XXX: currently this falls back to mpf_gamma. It should
# be the other way around: this function should handle
# all sizes/precisions efficiently, and gamma should fall
# back here
def mpf_gamma_int(n, prec, rounding=round_fast):
    if n < 1000:
        return from_int(int_fac(n-1), prec, rounding)
    # XXX: choose the cutoff less arbitrarily
    size = int(n*math.log(n,2))
    if prec > size/20.0:
        return from_int(int_fac(n-1), prec, rounding)
    return mpf_gamma(from_int(n), prec, rounding)

def mpf_factorial(x, prec, rounding=round_fast):
    return mpf_gamma(x, prec, rounding, p1=0)

def mpc_factorial(x, prec, rounding=round_fast):
    return mpc_gamma(x, prec, rounding, p1=0)

def mpf_gamma(x, prec, rounding=round_fast, p1=1):
    """
    Computes the gamma function of a real floating-point argument.
    With p1=0, computes a factorial instead.
    """
    sign, man, exp, bc = x
    if not man:
        if x == finf:
            return finf
        if x == fninf or x == fnan:
            return fnan
    # More precision is needed for enormous x. TODO:
    # use Stirling's formula + Euler-Maclaurin summation
    size = exp + bc
    if size > 5:
        size = int(size * math.log(size,2))
    wp = prec + max(0, size) + 15
    if exp >= 0:
        if sign or (p1 and not man):
            raise ValueError("gamma function pole")
        # A direct factorial is fastest
        if exp + bc <= 10:
            return from_int(int_fac((man<<exp)-p1), prec, rounding)
    reflect = sign or exp+bc < -1
    if p1:
        # Should be done exactly!
        x = mpf_sub(x, fone)
    # x < 0.25
    if reflect:
        # gamma = pi / (sin(pi*x) * gamma(1-x))
        wp += 15
        pix = mpf_mul(x, mpf_pi(wp), wp)
        t = mpf_sin_pi(x, wp)
        g = mpf_gamma(mpf_sub(fone, x), wp)
        return mpf_div(pix, mpf_mul(t, g, wp), prec, rounding)
    sprec, a, c = get_spouge_coefficients(wp)
    s = spouge_sum_real(x, sprec, a, c)
    # gamma = exp(log(x+a)*(x+0.5) - xpa) * s
    xpa = mpf_add(x, from_int(a), wp)
    logxpa = mpf_log(xpa, wp)
    xph = mpf_add(x, fhalf, wp)
    t = mpf_sub(mpf_mul(logxpa, xph, wp), xpa, wp)
    t = mpf_mul(mpf_exp(t, wp), s, prec, rounding)
    return t

def mpc_gamma(x, prec, rounding=round_fast, p1=1):
    re, im = x
    if im == fzero:
        return mpf_gamma(re, prec, rounding, p1), fzero
    # More precision is needed for enormous x.
    sign, man, exp, bc = re
    isign, iman, iexp, ibc = im
    if re == fzero:
        size = iexp+ibc
    else:
        size = max(exp+bc, iexp+ibc)
    if size > 5:
        size = int(size * math.log(size,2))
    reflect = sign or (exp+bc < -1)
    wp = prec + max(0, size) + 25
    # Near x = 0 pole (TODO: other poles)
    if p1:
        if size < -prec-5:
            return mpc_add_mpf(mpc_div(mpc_one, x, 2*prec+10), \
                mpf_neg(mpf_euler(2*prec+10)), prec, rounding)
        elif size < -5:
            wp += (-2*size)
    if p1:
        # Should be done exactly!
        re_orig = re
        re = mpf_sub(re, fone, bc+abs(exp)+2)
        x = re, im
    if reflect:
        # Reflection formula
        wp += 15
        pi = mpf_pi(wp), fzero
        pix = mpc_mul(x, pi, wp)
        t = mpc_sin_pi(x, wp)
        u = mpc_sub(mpc_one, x, wp)
        g = mpc_gamma(u, wp)
        w = mpc_mul(t, g, wp)
        return mpc_div(pix, w, wp)
    # Extremely close to the real line?
    # XXX: reflection formula
    if iexp+ibc < -wp:
        a = mpf_gamma(re_orig, wp)
        b = mpf_psi0(re_orig, wp)
        gamma_diff = mpf_div(a, b, wp)
        return mpf_pos(a, prec, rounding), mpf_mul(gamma_diff, im, prec, rounding)
    sprec, a, c = get_spouge_coefficients(wp)
    s = spouge_sum_complex(re, im, sprec, a, c)
    # gamma = exp(log(x+a)*(x+0.5) - xpa) * s
    repa = mpf_add(re, from_int(a), wp)
    logxpa = mpc_log((repa, im), wp)
    reph = mpf_add(re, fhalf, wp)
    t = mpc_sub(mpc_mul(logxpa, (reph, im), wp), (repa, im), wp)
    t = mpc_mul(mpc_exp(t, wp), s, prec, rounding)
    return t


#-----------------------------------------------------------------------#
#                                                                       #
#                         Polygamma functions                           #
#                                                                       #
#-----------------------------------------------------------------------#

"""
For all polygamma (psi) functions, we use the Euler-Maclaurin summation
formula. It looks slightly different in the m = 0 and m > 0 cases.

For m = 0, we have
                                 oo
                                ___   B
       (0)                1    \       2 k    -2 k
    psi   (z)  ~ log z + --- -  )    ------  z
                         2 z   /___  (2 k)!
                               k = 1

Experiment shows that the minimum term of the asymptotic series
reaches 2^(-p) when Re(z) > 0.11*p. So we simply use the recurrence
for psi (equivalent, in fact, to summing to the first few terms
directly before applying E-M) to obtain z large enough.

Since, very crudely, log z ~= 1 for Re(z) > 1, we can use
fixed-point arithmetic  (if z is extremely large, log(z) itself
is a sufficient approximation, so we can stop there already).

For Re(z) << 0, we could use recurrence, but this is of course
inefficient for large negative z, so there we use the
reflection formula instead.

For m > 0, we have

                  N - 1
                   ___
  ~~~(m)       [  \          1    ]         1            1
  psi   (z)  ~ [   )     -------- ] +  ---------- +  -------- +
               [  /___        m+1 ]           m+1           m
                  k = 1  (z+k)    ]    2 (z+N)       m (z+N)

      oo
     ___    B
    \        2 k   (m+1) (m+2) ... (m+2k-1)
  +  )     ------  ------------------------
    /___   (2 k)!            m + 2 k
    k = 1               (z+N)

where ~~~ denotes the function rescaled by 1/((-1)^(m+1) m!).

Here again N is chosen to make z+N large enough for the minimum
term in the last series to become smaller than eps.

TODO: the current estimation of N for m > 0 is *very suboptimal*.

TODO: implement the reflection formula for m > 0, Re(z) << 0.
It is generally a combination of multiple cotangents. Need to
figure out a reasonably simple way to generate these formulas
on the fly.

TODO: maybe use exact algorithms to compute psi for integral
and certain rational arguments, as this can be much more
efficient. (On the other hand, the availability of these
special values provides a convenient way to test the general
algorithm.)
"""

# Harmonic numbers are just shifted digamma functions
# We should calculate these exactly when x is an integer
# and when doing so is faster.

def mpf_harmonic(x, prec, rnd):
    if x in (fzero, fnan, finf):
        return x
    a = mpf_psi0(mpf_add(fone, x, prec+5), prec)
    return mpf_add(a, mpf_euler(prec+5, rnd), prec, rnd)

def mpc_harmonic(z, prec, rnd):
    if z[1] == fzero:
        return (mpf_harmonic(z[0], prec, rnd), fzero)
    a = mpc_psi0(mpc_add_mpf(z, fone, prec+5), prec)
    return mpc_add_mpf(a, mpf_euler(prec+5, rnd), prec, rnd)

def mpf_psi0(x, prec, rnd=round_fast):
    """
    Computation of the digamma function (psi function of order 0)
    of a real argument.
    """
    sign, man, exp, bc = x
    wp = prec + 10
    if not man:
        if x == finf: return x
        if x == fninf or x == fnan: return fnan
    if x == fzero or (exp >= 0 and sign):
        raise ValueError("polygamma pole")
    # Reflection formula
    if sign and exp+bc > 3:
        c, s = mpf_cos_sin_pi(x, wp)
        q = mpf_mul(mpf_div(c, s, wp), mpf_pi(wp), wp)
        p = mpf_psi0(mpf_sub(fone, x, wp), wp)
        return mpf_sub(p, q, prec, rnd)
    # The logarithmic term is accurate enough
    if (not sign) and bc + exp > wp:
        return mpf_log(mpf_sub(x, fone, wp), prec, rnd)
    # Initial recurrence to obtain a large enough x
    m = to_int(x)
    n = int(0.11*wp) + 2
    s = MP_ZERO
    x = to_fixed(x, wp)
    one = MP_ONE << wp
    if m < n:
        for k in xrange(m, n):
            s -= (one << wp) // x
            x += one
    x -= one
    # Logarithmic term
    s += to_fixed(mpf_log(from_man_exp(x, -wp, wp), wp), wp)
    # Endpoint term in Euler-Maclaurin expansion
    s += (one << wp) // (2*x)
    # Euler-Maclaurin remainder sum
    x2 = (x*x) >> wp
    t = one
    prev = 0
    k = 1
    while 1:
        t = (t*x2) >> wp
        bsign, bman, bexp, bbc = mpf_bernoulli(2*k, wp)
        offset = (bexp + 2*wp)
        if offset >= 0: term = (bman << offset) // (t*(2*k))
        else:           term = (bman >> (-offset)) // (t*(2*k))
        if k & 1: s -= term
        else:     s += term
        if k > 2 and term >= prev:
            break
        prev = term
        k += 1
    return from_man_exp(s, -wp, wp, rnd)

def mpc_psi0(z, prec, rnd=round_fast):
    """
    Computation of the digamma function (psi function of order 0)
    of a complex argument.
    """
    re, im = z
    # Fall back to the real case
    if im == fzero:
        return (mpf_psi0(re, prec, rnd), fzero)
    wp = prec + 20
    sign, man, exp, bc = re
    # Reflection formula
    if sign and exp+bc > 3:
        c = mpc_cos_pi(z, wp)
        s = mpc_sin_pi(z, wp)
        q = mpc_mul_mpf(mpc_div(c, s, wp), mpf_pi(wp), wp)
        p = mpc_psi0(mpc_sub(mpc_one, z, wp), wp)
        return mpc_sub(p, q, prec, rnd)
    # Just the logarithmic term
    if (not sign) and bc + exp > wp:
        return mpc_log(mpc_sub(z, mpc_one, wp), prec, rnd)
    # Initial recurrence to obtain a large enough z
    w = to_int(re)
    n = int(0.11*wp) + 2
    s = mpc_zero
    if w < n:
        for k in xrange(w, n):
            s = mpc_sub(s, mpc_reciprocal(z, wp), wp)
            z = mpc_add_mpf(z, fone, wp)
    z = mpc_sub(z, mpc_one, wp)
    # Logarithmic and endpoint term
    s = mpc_add(s, mpc_log(z, wp), wp)
    s = mpc_add(s, mpc_div(mpc_half, z, wp), wp)
    # Euler-Maclaurin remainder sum
    z2 = mpc_square(z, wp)
    t = mpc_one
    prev = mpc_zero
    k = 1
    eps = mpf_shift(fone, -wp+2)
    while 1:
        t = mpc_mul(t, z2, wp)
        bern = mpf_bernoulli(2*k, wp)
        term = mpc_mpf_div(bern, mpc_mul_int(t, 2*k, wp), wp)
        s = mpc_sub(s, term, wp)
        szterm = mpc_abs(term, 10)
        if k > 2 and mpf_le(szterm, eps):
            break
        prev = term
        k += 1
    return s

# Currently unoptimized
def mpf_psi(m, x, prec, rnd=round_fast):
    """
    Computation of the polygamma function of arbitrary integer order
    m >= 0, for a real argument x.
    """
    if m == 0:
        return mpf_psi0(x, prec, rnd=round_fast)
    return mpc_psi(m, (x, fzero), prec, rnd)[0]

def mpc_psi(m, z, prec, rnd=round_fast):
    """
    Computation of the polygamma function of arbitrary integer order
    m >= 0, for a complex argument z.
    """
    if m == 0:
        return mpc_psi0(z, prec, rnd)
    re, im = z
    wp = prec + 20
    sign, man, exp, bc = re
    if not man:
        if re == finf and im == fzero:
            return (fzero, fzero)
        if re == fnan:
            return fnan
    # Recurrence
    w = to_int(re)
    n = int(0.4*wp + 4*m)
    s = mpc_zero
    if w < n:
        for k in xrange(w, n):
            t = mpc_pow_int(z, -m-1, wp)
            s = mpc_add(s, t, wp)
            z = mpc_add_mpf(z, fone, wp)
    zm = mpc_pow_int(z, -m, wp)
    z2 = mpc_pow_int(z, -2, wp)
    # 1/m*(z+N)^m
    integral_term = mpc_div_mpf(zm, from_int(m), wp)
    s = mpc_add(s, integral_term, wp)
    # 1/2*(z+N)^(-(m+1))
    s = mpc_add(s, mpc_mul_mpf(mpc_div(zm, z, wp), fhalf, wp), wp)
    a = m + 1
    b = 2
    k = 1
    # Important: we want to sum up to the *relative* error,
    # not the absolute error, because psi^(m)(z) might be tiny
    magn = mpc_abs(s, 10)
    magn = magn[2]+magn[3]
    eps = mpf_shift(fone, magn-wp+2)
    while 1:
        zm = mpc_mul(zm, z2, wp)
        bern = mpf_bernoulli(2*k, wp)
        scal = mpf_mul_int(bern, a, wp)
        scal = mpf_div(scal, from_int(b), wp)
        term = mpc_mul_mpf(zm, scal, wp)
        s = mpc_add(s, term, wp)
        szterm = mpc_abs(term, 10)
        if k > 2 and mpf_le(szterm, eps):
            break
        #print k, to_str(szterm, 10), to_str(eps, 10)
        a *= (m+2*k)*(m+2*k+1)
        b *= (2*k+1)*(2*k+2)
        k += 1
    # Scale and sign factor
    v = mpc_mul_mpf(s, mpf_gamma(from_int(m+1), wp), prec, rnd)
    if not (m & 1):
        v = mpf_neg(v[0]), mpf_neg(v[1])
    return v


#-----------------------------------------------------------------------#
#                                                                       #
#                         Riemann zeta function                         #
#                                                                       #
#-----------------------------------------------------------------------#

"""
We use zeta(s) = eta(s) / (1 - 2**(1-s)) and Borwein's approximation

                  n-1
                  ___       k
             -1  \      (-1)  (d_k - d_n)
  eta(s) ~= ----  )     ------------------
             d_n /___              s
                 k = 0      (k + 1)
where
             k
             ___                i
            \     (n + i - 1)! 4
  d_k  =  n  )    ---------------.
            /___   (n - i)! (2i)!
            i = 0

If s = a + b*I, the absolute error for eta(s) is bounded by

    3 (1 + 2|b|)
    ------------ * exp(|b| pi/2)
               n
    (3+sqrt(8))

Disregarding the linear term, we have approximately,

  log(err) ~= log(exp(1.58*|b|)) - log(5.8**n)
  log(err) ~= 1.58*|b| - log(5.8)*n
  log(err) ~= 1.58*|b| - 1.76*n
  log2(err) ~= 2.28*|b| - 2.54*n

So for p bits, we should choose n > (p + 2.28*|b|) / 2.54.

References:
-----------

Peter Borwein, "An Efficient Algorithm for the Riemann Zeta Function"
http://www.cecm.sfu.ca/personal/pborwein/PAPERS/P117.ps

http://en.wikipedia.org/wiki/Dirichlet_eta_function
"""

borwein_cache = {}



def borwein_coefficients(n):
    if n in borwein_cache:
        return borwein_cache[n]
    ds = [MP_ZERO] * (n+1)
    d = MP_ONE
    s = ds[0] = MP_ONE
    for i in range(1, n+1):
        d = d * 4 * (n+i-1) * (n-i+1)
        d //= ((2*i) * ((2*i)-1))
        s += d
        ds[i] = s
    borwein_cache[n] = ds
    return ds

def mpf_zeta_int(s, prec, rnd=round_fast):
    """
    Optimized computation of zeta(s) for an integer s.
    """
    wp = prec + 20
    s = int(s)
    if s < 2:
        if s == 1:
            raise ValueError("zeta(1) pole")
        if not s:
            return mpf_neg(fhalf)
        return mpf_div(mpf_bernoulli(-s+1, wp), from_int(s-1), prec, rnd)
    # 2^-s term vanishes?
    if s >= wp:
        return mpf_perturb(fone, 0, prec, rnd)
    # 5^-s term vanishes?
    elif s >= wp*0.431:
        t = one = 1 << wp
        t += 1 << (wp - s)
        t += one // (MP_THREE ** s)
        t += 1 << max(0, wp - s*2)
        return from_man_exp(t, -wp, prec, rnd)
    else:
        # Fast enough to sum directly?
        # Even better, we use the Euler product (idea stolen from pari)
        m = (float(wp)/(s-1) + 1)
        if m < 30:
            needed_terms = int(2.0**m + 1)
            if needed_terms < int(wp/2.54 + 5) / 10:
                t = fone
                for k in list_primes(needed_terms):
                    #print k, needed_terms
                    powprec = int(wp - s*math.log(k,2))
                    if powprec < 2:
                        break
                    a = mpf_sub(fone, mpf_pow_int(from_int(k), -s, powprec), wp)
                    t = mpf_mul(t, a, wp)
                return mpf_div(fone, t, wp)
    # Use Borwein's algorithm
    n = int(wp/2.54 + 5)
    d = borwein_coefficients(n)
    t = MP_ZERO
    s = MP_BASE(s)
    for k in xrange(n):
        t += (((-1)**k * (d[k] - d[n])) << wp) // (k+1)**s
    t = (t << wp) // (-d[n])
    t = (t << wp) // ((1 << wp) - (1 << (wp+1-s)))
    return from_man_exp(t, -wp-wp, prec, rnd)

def mpf_zeta(s, prec, rnd=round_fast, alt=0):
    sign, man, exp, bc = s
    if not man:
        if s == fzero:
            if alt:
                return fhalf
            else:
                return mpf_neg(fhalf)
        if s == finf:
            return fone
        return fnan
    wp = prec + 20
    # First term vanishes?
    if (not sign) and (exp + bc > (math.log(wp,2) + 2)):
        return mpf_perturb(fone, alt, prec, rnd)
    # Optimize for integer arguments
    elif exp >= 0:
        if alt:
            if s == fone:
                return mpf_ln2(prec, rnd)
            z = mpf_zeta_int(to_int(s), wp, negative_rnd[rnd])
            q = mpf_sub(fone, mpf_pow(ftwo, mpf_sub(fone, s, wp), wp), wp)
            return mpf_mul(z, q, prec, rnd)
        else:
            return mpf_zeta_int(to_int(s), prec, rnd)
    # Negative: use the reflection formula
    # Borwein only proves the accuracy bound for x >= 1/2. However, based on
    # tests, the accuracy without reflection is quite good even some distance
    # to the left of 1/2. XXX: verify this.
    if sign:
        # XXX: could use the separate refl. formula for Dirichlet eta
        if alt:
            q = mpf_sub(fone, mpf_pow(ftwo, mpf_sub(fone, s, wp), wp), wp)
            return mpf_mul(mpf_zeta(s, wp), q, prec, rnd)
        # XXX: -1 should be done exactly
        y = mpf_sub(fone, s, 10*wp)
        a = mpf_gamma(y, wp)
        b = mpf_zeta(y, wp)
        c = mpf_sin_pi(mpf_shift(s, -1), wp)
        wp2 = wp + (exp+bc)
        pi = mpf_pi(wp+wp2)
        d = mpf_div(mpf_pow(mpf_shift(pi, 1), s, wp2), pi, wp2)
        return mpf_mul(a,mpf_mul(b,mpf_mul(c,d,wp),wp),prec,rnd)
    t = MP_ZERO
    #wp += 16 - (prec & 15)
    # Use Borwein's algorithm
    n = int(wp/2.54 + 5)
    d = borwein_coefficients(n)
    t = MP_ZERO
    sf = to_fixed(s, wp)
    for k in xrange(n):
        u = from_man_exp(-sf*log_int_fixed(k+1, wp), -2*wp, wp)
        esign, eman, eexp, ebc = mpf_exp(u, wp)
        offset = eexp + wp
        if offset >= 0:
            w = ((d[k] - d[n]) * eman) << offset
        else:
            w = ((d[k] - d[n]) * eman) >> (-offset)
        if k & 1:
            t -= w
        else:
            t += w
    t = t // (-d[n])
    t = from_man_exp(t, -wp, wp)
    if alt:
        return mpf_pos(t, prec, rnd)
    else:
        q = mpf_sub(fone, mpf_pow(ftwo, mpf_sub(fone, s, wp), wp), wp)
        return mpf_div(t, q, prec, rnd)

def mpc_zeta(s, prec, rnd=round_fast, alt=0):
    re, im = s
    if im == fzero:
        return mpf_zeta(re, prec, rnd, alt), fzero
    wp = prec + 20
    # Reflection formula. To be rigorous, we should reflect to the left of
    # re = 1/2 (see comments for mpf_zeta), but this leads to unnecessary
    # slowdown for interesting values of s
    if mpf_lt(re, fzero):
        # XXX: could use the separate refl. formula for Dirichlet eta
        if alt:
            q = mpc_sub(mpc_one, mpc_pow(mpc_two, mpc_sub(mpc_one, s, wp),
                wp), wp)
            return mpc_mul(mpc_zeta(s, wp), q, prec, rnd)
        # XXX: -1 should be done exactly
        y = mpc_sub(mpc_one, s, 10*wp)
        a = mpc_gamma(y, wp)
        b = mpc_zeta(y, wp)
        c = mpc_sin_pi(mpc_shift(s, -1), wp)
        rsign, rman, rexp, rbc = re
        isign, iman, iexp, ibc = im
        mag = max(rexp+rbc, iexp+ibc)
        wp2 = wp + mag
        pi = mpf_pi(wp+wp2)
        pi2 = (mpf_shift(pi, 1), fzero)
        d = mpc_div_mpf(mpc_pow(pi2, s, wp2), pi, wp2)
        return mpc_mul(a,mpc_mul(b,mpc_mul(c,d,wp),wp),prec,rnd)
    n = int(wp/2.54 + 5)
    n += int(0.9*abs(to_int(im)))
    d = borwein_coefficients(n)
    ref = to_fixed(re, wp)
    imf = to_fixed(im, wp)
    tre = MP_ZERO
    tim = MP_ZERO
    one = MP_ONE << wp
    one_2wp = MP_ONE << (2*wp)
    critical_line = re == fhalf
    for k in xrange(n):
        log = log_int_fixed(k+1, wp)
        # A square root is much cheaper than an exp
        if critical_line:
            w = one_2wp // sqrt_fixed((k+1) << wp, wp)
        else:
            w = to_fixed(mpf_exp(from_man_exp(-ref*log, -2*wp), wp), wp)
        if k & 1:
            w *= (d[n] - d[k])
        else:
            w *= (d[k] - d[n])
        wre, wim = cos_sin(from_man_exp(-imf * log_int_fixed(k+1, wp), -2*wp), wp)
        tre += (w * to_fixed(wre, wp)) >> wp
        tim += (w * to_fixed(wim, wp)) >> wp
    tre //= (-d[n])
    tim //= (-d[n])
    tre = from_man_exp(tre, -wp, wp)
    tim = from_man_exp(tim, -wp, wp)
    if alt:
        return mpc_pos((tre, tim), prec, rnd)
    else:
        q = mpc_sub(mpc_one, mpc_pow(mpc_two, mpc_sub(mpc_one, s, wp), wp), wp)
        return mpc_div((tre, tim), q, prec, rnd)

def mpf_altzeta(s, prec, rnd=round_fast):
    return mpf_zeta(s, prec, rnd, 1)

def mpc_altzeta(s, prec, rnd=round_fast):
    return mpc_zeta(s, prec, rnd, 1)
