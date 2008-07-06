"""
Miscellaneous special functions
"""

from lib import *
from libmpc import *
from mptypes import *
from mptypes import constant

__docformat__ = 'plaintext'

#---------------------------------------------------------------------------#
#                                                                           #
#                       First some mathematical constants                   #
#                                                                           #
#---------------------------------------------------------------------------#


# The golden ratio is given by phi = (1 + sqrt(5))/2

@constant_memo
def phi_fixed(prec):
    prec += 10
    sqrt = [sqrt_fixed2, sqrt_fixed][prec < 20000]
    a = sqrt(MP_FIVE<<prec, prec) + (MP_ONE << prec)
    return a >> 11

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

# Euler's constant (gamma) is computed using the Brent-McMillan formula,
# gamma ~= A(n)/B(n) - log(n), where

#   A(n) = sum_{k=0,1,2,...} (n**k / k!)**2 * H(k)
#   B(n) = sum_{k=0,1,2,...} (n**k / k!)**2
#   H(k) = 1 + 1/2 + 1/3 + ... + 1/k

# The error is bounded by O(exp(-4n)). Choosing n to be a power
# of two, 2**p, the logarithm becomes particularly easy to calculate.

# Reference:
# Xavier Gourdon & Pascal Sebah, The Euler constant: gamma
# http://numbers.computation.free.fr/Constants/Gamma/gamma.pdf

@constant_memo
def euler_fixed(prec):
    prec += 30
    # choose p such that exp(-4*(2**p)) < 2**-n
    p = int(math.log((prec/4) * math.log(2), 2)) + 1
    n = MP_ONE<<p
    r = one = MP_ONE<<prec
    H, A, B, npow, k, d = MP_ZERO, MP_ZERO, MP_ZERO, 1, 1, 1
    while r:
        A += (r * H) >> prec
        B += r
        r = r * (n*n) // (k*k)
        H += one // k
        k += 1
    S = ((A<<prec) // B) - p*log2_fixed(prec)
    return S >> 30

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
    orig = mp.prec
    try:
        mp.prec = int(prec + prec**0.5 + 15)
        s = mpf(0)
        t = one = mpf(1)
        B = bernoulli_range()
        fac = mpf(4)
        pipow = twopi2 = (2*pi)**2
        n = 1
        while 1:
            zeta2n = (-1)**(n+1) * B.next() * pipow / fac
            term = ((zeta2n - 1) * t) / n
            # print n, nstr(term)
            if term < eps:
                break
            s += term
            t += (one/(2*n+1) - one/(2*n))
            n += 1
            fac *= (2*n)*(2*n-1)
            pipow *= twopi2
        return to_fixed(exp(s/ln2)._mpf_, prec)
    finally:
        mp.prec = orig

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
    orig = mp.prec
    try:
        dps = mp.dps
        mp.prec = prec + 30
        N = int(1.0*dps + 5)
        logs = log_range()
        s = mpf(0)
        # E-M step 1: sum log(k)/k**2 for k = 2..N-1
        for n in range(2, N):
            # print n, N
            logn = logs.next()
            s += logn / n**2
        logN = logs.next()
        # E-M step 2: integral of log(x)/x**2 from N to inf
        s += (1+logN)/N
        # E-M step 3: endpoint correction term f(N)/2
        s += logN/(N**2 * 2)
        # E-M step 4: the series of derivatives
        pN, a, b, j, B2k, fac, k = N**3, 1, -2, 3, bernoulli_range(), 2, 1
        while 1:
            # D(2*k-1) * B(2*k) / fac(2*k) [D(n) = nth derivative]
            D = (a+b*logN)/pN
            B = B2k.next()
            term = B * D / fac
            if abs(term) < eps:
                break
            # print k, nstr(term)
            s -= term
            # Advance derivative twice
            a, b, pN, j = b-a*j, -j*b, pN*N, j+1
            a, b, pN, j = b-a*j, -j*b, pN*N, j+1
            k += 1
            fac *= (2*k) * (2*k-1)
        A = exp((6*s/pi**2 + log(2*pi) + euler)/12)
        return to_fixed(A._mpf_, prec)
    finally:
        mp.prec = orig

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

fme = from_man_exp

def mpf_phi(p, r): return fme(phi_fixed(p+10), -p-10, p, r)
def mpf_khinchin(p, r): return fme(khinchin_fixed(p+10), -p-10, p, r)
def mpf_glaisher(p, r): return fme(glaisher_fixed(p+10), -p-10, p, r)
def mpf_apery(p, r): return fme(apery_fixed(p+10), -p-10, p, r)
def mpf_euler(p, r): return fme(euler_fixed(p+10), -p-10, p, r)
def mpf_catalan(p, r): return fme(catalan_fixed(p+10), -p-10, p, r)

phi = constant(mpf_phi, "Golden ratio (phi)")
catalan = constant(mpf_catalan, "Catalan's constant")
euler = constant(mpf_euler, "Euler's constant (gamma)")
khinchin = constant(mpf_khinchin, "Khinchin's constant")
glaisher = constant(mpf_glaisher, "Glaisher's constant")
apery = constant(mpf_apery, "Apery's constant")


#----------------------------------------------------------------------
# Factorial related functions
#

# For internal use
def int_fac(n, memo={0:1, 1:1}):
    """Return n factorial (for integers n >= 0 only)."""
    f = memo.get(n)
    if f:
        return f
    k = len(memo)
    p = memo[k-1]
    while k <= n:
        p *= k
        if k < 1024:
            memo[k] = p
        k += 1
    return p

if MODE == "gmpy":
    int_fac = gmpy.fac

"""
We compute the gamma function using Spouge's approximation

    x! = (x+a)**(x+1/2) * exp(-x-a) * [c_0 + S(x) + eps]

where S(x) is the sum of c_k/(x+k) from k = 1 to a-1 and the coefficients
are given by

  c_0 = sqrt(2*pi)

         (-1)**(k-1)
  c_k =  ----------- (a-k)**(k-1/2) exp(-k+a),  k = 1,2,...,a-1
          (k - 1)!

Due to an inequality proved by Spouge, if we choose a = int(1.26*n), the
error eps is less than 10**-n for any x in the right complex half-plane
(assuming a > 2). In practice, it seems that a can be chosen quite a bit
lower still (30-50%); this possibility should be investigated.

Reference:
John L. Spouge, "Computation of the gamma, digamma, and trigamma
functions", SIAM Journal on Numerical Analysis 31 (1994), no. 3, 931-944.
"""

spouge_cache = {}

def calc_spouge_coefficients(a, prec):
    wp = prec + int(a*1.4)
    c = [0] * a
    # b = exp(a-1)
    b = fexp(from_int(a-1), wp)
    # e = exp(1)
    e = fexp(fone, wp)
    # sqrt(2*pi)
    sq2pi = fsqrt(fshift(fpi(wp), 1), wp)
    c[0] = to_fixed(sq2pi, prec)
    for k in xrange(1, a):
        # c[k] = ((-1)**(k-1) * (a-k)**k) * b / sqrt(a-k)
        term = fmuli(b, ((-1)**(k-1) * (a-k)**k), wp)
        term = fdiv(term, fsqrt(from_int(a-k), wp), wp)
        c[k] = to_fixed(term, prec)
        # b = b / (e * k)
        b = fdiv(b, fmul(e, from_int(k), wp), wp)
    return c

# Cached lookup of coefficients
def get_spouge_coefficients(prec):

    # This exact precision has been used before
    if prec in spouge_cache:
        return spouge_cache[prec]

    for p in spouge_cache:
        if 0.8 <= float(p)/prec < 1:
            return spouge_cache[p]

    # Here we estimate the value of a based on Spouge's inequality for
    # the relative error
    a = max(3, int(0.39*prec))  # ~= 1.26*n

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
#
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

def mpf_gamma(x, prec, rounding=round_fast, p1=1):
    sign, man, exp, bc = x
    if exp >= 0:
        if sign or (p1 and not man):
            raise ValueError("gamma function pole")
        # A direct factorial is fastest
        if exp + bc <= 10:
            return from_int(int_fac((man<<exp)-p1), prec, rounding)
    wp = prec + 15
    if p1:
        x = fsub(x, fone, wp)
    # x < 0.25
    if sign or exp+bc < -1:
        # gamma = pi / (sin(pi*x) * gamma(1-x))
        wp += 15
        pi = fpi(wp)
        pix = fmul(x, pi, wp)
        t = fsin(pix, wp)
        g = mpf_gamma(fsub(fone, x, wp), wp)
        return fdiv(pix, fmul(t, g, wp), prec, rounding)
    sprec, a, c = get_spouge_coefficients(wp)
    s = spouge_sum_real(x, sprec, a, c)
    # gamma = exp(log(x+a)*(x+0.5) - xpa) * s
    xpa = fadd(x, from_int(a), wp)
    logxpa = flog(xpa, wp)
    xph = fadd(x, fhalf, wp)
    t = fsub(fmul(logxpa, xph, wp), xpa, wp)
    t = fmul(fexp(t, wp), s, prec, rounding)
    return t

def mpc_gamma(x, prec, rounding=round_fast, p1=1):
    re, im = x
    if im == fzero:
        return mpf_gamma(re, prec, rounding, p1), fzero
    wp = prec + 25
    sign, man, exp, bc = re
    if p1:
        re = fsub(re, fone, wp)
        x = re, im
    if sign or exp+bc < -1:
        # Reflection formula
        wp += 15
        pi = fpi(wp), fzero
        pix = mpc_mul(x, pi, wp)
        t = mpc_sin(pix, wp)
        u = mpc_sub(mpc_one, x, wp)
        g = mpc_gamma(u, wp)
        w = mpc_mul(t, g, wp)
        return mpc_div(pix, w, wp)
    sprec, a, c = get_spouge_coefficients(wp)
    s = spouge_sum_complex(re, im, sprec, a, c)
    # gamma = exp(log(x+a)*(x+0.5) - xpa) * s
    repa = fadd(re, from_int(a), wp)
    logxpa = mpc_log((repa, im), wp)
    reph = fadd(re, fhalf, wp)
    t = mpc_sub(mpc_mul(logxpa, (reph, im), wp), (repa, im), wp)
    t = mpc_mul(mpc_exp(t, wp), s, prec, rounding)
    return t

def gamma(x):
    x = convert_lossless(x)
    prec = mp.prec
    if isinstance(x, mpf):
        return make_mpf(mpf_gamma(x._mpf_, prec, round_nearest, 1))
    else:
        return make_mpc(mpc_gamma(x._mpc_, prec, round_nearest, 1))

def factorial(x):
    x = convert_lossless(x)
    prec = mp.prec
    if isinstance(x, mpf):
        return make_mpf(mpf_gamma(x._mpf_, prec, round_nearest, 0))
    else:
        return make_mpc(mpc_gamma(x._mpc_, prec, round_nearest, 0))

def isnpint(x):
    if not x:
        return True
    if isinstance(x, mpf):
        sign, man, exp, bc = x._mpf_
        return sign and exp >= 0
    if isinstance(x, mpc):
        return not x.imag and isnpint(x.real)

def gammaprod(a, b):
    """
    Computes the product / quotient of gamma functions

        G(a_0) G(a_1) ... G(a_p)
        ------------------------
        G(b_0) G(b_1) ... G(a_q)

    with proper cancellation of poles (interpreting the expression as a
    limit). Returns +inf if the limit diverges.
    """
    a = [convert_lossless(x) for x in a]
    b = [convert_lossless(x) for x in b]
    poles_num = []
    poles_den = []
    regular_num = []
    regular_den = []
    for x in a: [regular_num, poles_num][isnpint(x)].append(x)
    for x in b: [regular_den, poles_den][isnpint(x)].append(x)
    # One more pole in numerator or denominator gives 0 or inf
    if len(poles_num) < len(poles_den): return mpf(0)
    if len(poles_num) > len(poles_den): return mpf('+inf')
    # All poles cancel
    # lim G(i)/G(j) = (-1)**(i+j) * gamma(1-j) / gamma(1-i)
    p = mpf(1)
    orig = mp.prec
    try:
        mp.prec = orig + 15
        while poles_num:
            i = poles_num.pop()
            j = poles_den.pop()
            p *= (-1)**(i+j) * gamma(1-j) / gamma(1-i)
        for x in regular_num: p *= gamma(x)
        for x in regular_den: p /= gamma(x)
    finally:
        mp.prec = orig
    return +p

def binomial(n, k):
    """Binomial coefficient, C(n,k) = n!/(k!*(n-k)!)."""
    return gammaprod([n+1], [k+1, n-k+1])

def rf(x, n):
    """Rising factorial (Pochhammer symbol), x^(n)"""
    return gammaprod([x+n], [x])

def ff(x, n):
    """Falling factorial, x_(n)"""
    return gammaprod([x+1], [x-n+1])



#---------------------------------------------------------------------------#
#                                                                           #
#                         Riemann zeta function                             #
#                                                                           #
#---------------------------------------------------------------------------#

"""
We use zeta(s) = eta(s) * (1 - 2**(1-s)) and Borwein's approximation
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

Reference:
Peter Borwein, "An Efficient Algorithm for the Riemann Zeta Function"
http://www.cecm.sfu.ca/personal/pborwein/PAPERS/P117.ps

http://en.wikipedia.org/wiki/Dirichlet_eta_function
"""

d_cache = {}

def zeta_coefs(n):
    if n in d_cache:
        return d_cache[n]
    ds = [MP_ZERO] * (n+1)
    d = MP_ONE
    s = ds[0] = MP_ONE
    for i in range(1, n+1):
        d = d * 4 * (n+i-1) * (n-i+1)
        d //= ((2*i) * ((2*i)-1))
        s += d
        ds[i] = s
    d_cache[n] = ds
    return ds

# Integer logarithms
_log_cache = {}

def _logk(k):
    p = mp.prec
    if k in _log_cache and _log_cache[k][0] >= p:
        return +_log_cache[k][1]
    else:
        x = log(k)
        _log_cache[k] = (p, x)
        return x

@extraprec(10, normalize_output=True)
def zeta(s):
    """Returns the Riemann zeta function of s."""
    s = convert_lossless(s)
    if s.real < 0:
        # Reflection formula (XXX: gets bad around the zeros)
        return 2**s * pi**(s-1) * sin(pi*s/2) * gamma(1-s) * zeta(1-s)
    else:
        p = mp.prec
        n = int((p + 2.28*abs(float(mpc(s).imag)))/2.54) + 3
        d = zeta_coefs(n)
        if isinstance(s, mpf) and s == int(s):
            sint = int(s)
            t = 0
            for k in range(n):
                t += (((-1)**k * (d[k] - d[n])) << p) // (k+1)**sint
            return (mpf((t, -p)) / -d[n]) / (1 - mpf(2)**(1-sint))
        else:
            t = mpf(0)
            for k in range(n):
                t += (-1)**k * mpf(d[k]-d[n]) * exp(-_logk(k+1)*s)
            return (t / -d[n]) / (mpf(1) - exp(log(2)*(1-s)))


@extraprec(5, normalize_output=True)
def bernoulli(n):
    """nth Bernoulli number, B_n"""
    if n == 1:
        return mpf(-0.5)
    if n & 1:
        return mpf(0)
    m = n // 2
    return (-1)**(m-1) * 2 * factorial(n) / (2*pi)**n * zeta(n)

# For sequential computation of Bernoulli numbers, we use Ramanujan's formula

#                            / n + 3 \
#   B   =  (A(n) - S(n))  /  |       |
#    n                       \   n   /

# where A(n) = (n+3)/3 when n = 0 or 2 (mod 6), A(n) = -(n+3)/6
# when n = 4 (mod 6), and

#          [n/6]
#           ___
#          \      /  n + 3  \
#   S(n) =  )     |         | * B
#          /___   \ n - 6*k /    n-6*k
#          k = 1

def bernoulli_range():
    """Generates B(2), B(4), B(6), ..."""
    oprec = mp.prec
    rounding = mp.rounding[0]
    prec = oprec + 30
    computed = {0:fone}
    m, bin1, bin = 2, MP_ONE, MP_BASE(10)
    f3 = from_int(3)
    f6 = from_int(6)
    while 1:
        case = m % 6
        s = fzero
        if m < 6: a = MP_ZERO
        else:     a = bin1
        for j in xrange(1, m//6+1):
            s = fadd(s, fmuli(computed[m-6*j], a, prec), prec)
            # Inner binomial coefficient
            j6 = 6*j
            a *= ((m-5-j6)*(m-4-j6)*(m-3-j6)*(m-2-j6)*(m-1-j6)*(m-j6))
            a //= ((4+j6)*(5+j6)*(6+j6)*(7+j6)*(8+j6)*(9+j6))
        if case == 0: b = fdivi(m+3, f3, prec)
        if case == 2: b = fdivi(m+3, f3, prec)
        if case == 4: b = fdivi(-m-3, f6, prec)
        b = fdiv(fsub(b, s, prec), from_int(bin), prec)
        computed[m] = b
        yield make_mpf(fpos(b, oprec, rounding))
        m += 2
        bin = bin * ((m+2)*(m+3)) // (m*(m-1))
        if m > 6: bin1 = bin1 * ((2+m)*(3+m)) // ((m-7)*(m-6))


#---------------------------------------------------------------------------#
#                                                                           #
#                          Hypergeometric functions                         #
#                                                                           #
#---------------------------------------------------------------------------#

import operator

"""
TODO:
  * By counting the number of multiplications vs divisions,
    the bit size of p can be kept around wp instead of growing
    it to n*wp for some (possibly large) n

  * Due to roundoff error, the series may fail to converge
    when x is negative and the convergence is slow.

"""

def hypsum(ar, af, ac, br, bf, bc, x):
    """
    Generic hypergeometric summation. This function computes:

            1   a_1 a_2 ...     1  (a_1 + 1) (a_2 + 1) ...  2
        1 + --  ----------- x + -- ----------------------- x  + ...
            1!  b_1 b_2 ...     2! (b_1 + 1) (b_2 + 1) ...

    The a_i and b_i sequences are separated by type:

    ar - list of a_i rationals [p,q]
    af - list of a_i mpf value tuples
    ac - list of a_i mpc value tuples
    br - list of b_i rationals [p,q]
    bf - list of b_i mpf value tuples
    bc - list of b_i mpc value tuples

    Note: the rational coefficients will be updated in-place and must
    hence be mutable (lists rather than tuples).

    x must be an mpf or mpc instance.
    """

    have_float = af or bf
    have_complex = ac or bc

    prec = mp.prec
    rnd = mp.rounding[0]
    wp = prec + 25

    if isinstance(x, mpf):
        x = to_fixed(x._mpf_, wp)
        y = MP_ZERO
    else:
        have_complex = 1
        x, y = x._mpc_
        x = to_fixed(x, wp)
        y = to_fixed(y, wp)

    sre = pre = one = MP_ONE << wp
    sim = pim = MP_ZERO

    n = 1

    # Need to shift down by wp once for each fixed-point multiply
    # At minimum, we multiply by once by x each step
    shift = 1

    # Fixed-point real coefficients
    if have_float:
        len_af = len(af)
        len_bf = len(bf)
        range_af = range(len_af)
        range_bf = range(len_bf)
        for i in range_af: af[i] = to_fixed(af[i], wp)
        for i in range_bf: bf[i] = to_fixed(bf[i], wp)
        shift += len_af

    if have_complex:
        len_ac = len(ac)
        len_bc = len(bc)
        range_ac = range(len_ac)
        range_bc = range(len_bc)
        for i in range_ac: ac[i] = [to_fixed(ac[i][0], wp), to_fixed(ac[i][1], wp)]
        for i in range_bc: bc[i] = [to_fixed(bc[i][0], wp), to_fixed(bc[i][1], wp)]
        shift += len_ac

    aqs = [a[1] for a in ar]
    bqs = [b[1] for b in br]
    aqprod = reduce(operator.mul, aqs, 1)
    bqprod = reduce(operator.mul, bqs, 1)

    assert shift >= 0

    while 1:
        # Integer and rational part of product
        mul = bqprod
        div = n * aqprod
        for ap, aq in ar: mul *= ap
        for bp, bq in br: div *= bp

        if have_complex:
            # Multiply by rational factors
            pre *= mul
            pim *= mul
            # Multiply by z
            pre, pim = pre*x - pim*y, pim*x + pre*y
            # Multiply by real factors
            for a in af:
                pre *= a
                pim *= a
            # Multiply by complex factors
            for are, aim in ac:
                pre, pim = pre*are - pim*aim, pim*are + pre*aim
            # Divide by rational factors
            pre //= div
            pim //= div
            # Divide by real factors
            for b in bf:
                pre = (pre << wp) // b
                pim = (pim << wp) // b
            # Divide by complex factors
            for bre, bim in bc:
                mag = bre*bre + bim*bim
                re = pre*bre + pim*bim
                im = pim*bre - pre*bim
                pre = (re << wp) // mag
                pim = (im << wp) // mag
        elif have_float:
            # Multiply and divide by real and rational factors, and x
            for a in af: pre *= a
            for b in bf:
                pre = (pre << wp) // b
            pre = (pre * (mul * x)) // div

        else:
            # Multiply and divide by rational factors and x
            pre = (pre * (mul * x)) // div

        pre >>= (wp*shift)
        sre += pre

        if have_complex:
            pim >>= (wp*shift)
            sim += pim
            if (-100 < pre < 100) and (-100 < pim < 100):
                break
        else:
            if -100 < pre < 100:
                break

        # Add 1 to all as and bs
        n += 1
        for ap_aq in ar: ap_aq[0] += ap_aq[1]
        for bp_bq in br: bp_bq[0] += bp_bq[1]
        if have_float:
            for i in range_af: af[i] += one
            for i in range_bf: bf[i] += one
        if have_complex:
            for i in range_ac: ac[i][0] += one
            for i in range_bc: bc[i][0] += one

    re = from_man_exp(sre, -wp, prec, rnd)
    if have_complex:
        return make_mpc((re, from_man_exp(sim, -wp, prec, rnd)))
    else:
        return make_mpf(re)


#---------------------------------------------------------------------------#
#   Special-case implementation for rational parameters. These are          #
#   about 2x faster at low precision                                        #
#---------------------------------------------------------------------------#

def sum_hyp0f1_rat((bp, bq), x):
    """Sum 0F1 for rational a. x must be mpf or mpc."""
    prec = mp.prec
    rnd = mp.rounding[0]
    wp = prec + 25
    if isinstance(x, mpf):
        x = to_fixed(x._mpf_, wp)
        s = p = MP_ONE << wp
        n = 1
        while 1:
            p = (p * (bq*x) // (n*bp)) >> wp
            if -100 < p < 100:
                break
            s += p; n += 1; bp += bq
        return make_mpf(from_man_exp(s, -wp, prec, rnd))
    else:
        wp = prec + 25
        zre, zim = x._mpc_
        zre = to_fixed(zre, wp)
        zim = to_fixed(zim, wp)
        sre = pre = MP_ONE << wp
        sim = pim = MP_ZERO
        n = 1
        while 1:
            r1 = bq
            r2 = n*bp
            pre, pim = pre*zre - pim*zim, pim*zre + pre*zim
            pre = ((pre * r1) // r2) >> wp
            pim = ((pim * r1) // r2) >> wp
            if -100 < pre < 100 and -100 < pim < 100:
                break
            sre += pre; sim += pim; n += 1; bp += bq
        re = from_man_exp(sre, -wp, prec, rnd)
        im = from_man_exp(sim, -wp, prec, rnd)
        return make_mpc((re, im))


def sum_hyp1f1_rat((ap, aq), (bp, bq), x):
    """Sum 1F1 for rational a, b. x must be mpf or mpc."""
    prec = mp.prec
    rnd = mp.rounding[0]
    wp = prec + 25
    if isinstance(x, mpf):
        x = to_fixed(x._mpf_, wp)
        s = p = MP_ONE << wp
        n = 1
        while 1:
            p = (p * (ap*bq*x) // (n*aq*bp)) >> wp
            if -100 < p < 100:
                break
            s += p; n += 1; ap += aq; bp += bq
        return make_mpf(from_man_exp(s, -wp, prec, rnd))
    else:
        wp = prec + 25
        zre, zim = x._mpc_
        zre = to_fixed(zre, wp)
        zim = to_fixed(zim, wp)
        sre = pre = MP_ONE << wp
        sim = pim = MP_ZERO
        n = 1
        while 1:
            r1 = ap*bq
            r2 = n*aq*bp
            pre, pim = pre*zre - pim*zim, pim*zre + pre*zim
            pre = ((pre * r1) // r2) >> wp
            pim = ((pim * r1) // r2) >> wp
            if -100 < pre < 100 and -100 < pim < 100:
                break
            sre += pre; sim += pim; n += 1; ap += aq; bp += bq
        re = from_man_exp(sre, -wp, prec, rnd)
        im = from_man_exp(sim, -wp, prec, rnd)
        return make_mpc((re, im))

def sum_hyp2f1_rat((ap, aq), (bp, bq), (cp, cq), x):
    """Sum 2F1 for rational a, b, c. x must be mpf or mpc"""
    prec = mp.prec
    rnd = mp.rounding[0]
    wp = prec + 25
    if isinstance(x, mpf):
        x = to_fixed(x._mpf_, wp)
        s = p = MP_ONE << wp
        n = 1
        while 1:
            p = (p * (ap*bp*cq*x) // (n*aq*bq*cp)) >> wp
            if -100 < p < 100:
                break
            s += p; n += 1; ap += aq; bp += bq; cp += cq
        return make_mpf(from_man_exp(s, -wp, prec, rnd))
    else:
        wp = prec + 25
        zre, zim = x._mpc_
        zre = to_fixed(zre, wp)
        zim = to_fixed(zim, wp)
        sre = pre = MP_ONE << wp
        sim = pim = MP_ZERO
        n = 1
        while 1:
            r1 = ap*bp*cq
            r2 = n*aq*bq*cp
            pre, pim = pre*zre - pim*zim, pim*zre + pre*zim
            pre = ((pre * r1) // r2) >> wp
            pim = ((pim * r1) // r2) >> wp
            if -100 < pre < 100 and -100 < pim < 100:
                break
            sre += pre; sim += pim; n += 1; ap += aq; bp += bq; cp += cq
        re = from_man_exp(sre, -wp, prec, rnd)
        im = from_man_exp(sim, -wp, prec, rnd)
        return make_mpc((re, im))

def parse_param(x):
    if isinstance(x, tuple):
        p, q = x
        return [[p, q]], [], []
    if isinstance(x, (int, long)):
        return [[x, 1]], [], []
    x = convert_lossless(x)
    if isinstance(x, mpf):
        return [], [x._mpf_], []
    if isinstance(x, mpc):
        return [], [], [x._mpc_]

class _mpq(tuple):
    @property
    def _mpf_(self):
        return (mpf(self[0])/self[1])._mpf_
    def __add__(self, other):
        if isinstance(other, _mpq):
            a, b = self
            c, d = other
            return _mpq((a*d+b*c, b*d))
        return NotImplemented
    def __sub__(self, other):
        if isinstance(other, _mpq):
            a, b = self
            c, d = other
            return _mpq((a*d-b*c, b*d))
        return NotImplemented

_1 = _mpq((1,1))
_0 = _mpq((0,1))

def _as_num(x):
    if isinstance(x, list):
        return _mpq(x)
    return x

def eval_hyp2f1(a,b,c,z):
    ar, af, ac = parse_param(a)
    br, bf, bc = parse_param(b)
    cr, cf, cc = parse_param(c)
    absz = abs(z)
    if absz == 1:
        # TODO: determine whether it actually does, and otherwise
        # return infinity instead
        print "Warning: 2F1 might not converge for |z| = 1"
    if absz <= 1:
        if ar and br and cr:
            return sum_hyp2f1_rat(ar[0], br[0], cr[0], z)
        return hypsum(ar+br, af+bf, ac+bc, cr, cf, cc, z)
    # Use 1/z transformation
    a = (ar and _as_num(ar[0])) or convert_lossless(a)
    b = (br and _as_num(br[0])) or convert_lossless(b)
    c = (cr and _as_num(cr[0])) or convert_lossless(c)
    orig = mp.prec
    try:
        mp.prec = orig + 15
        h1 = eval_hyp2f1(a, _1-c+a, _1-b+a, 1/z)
        h2 = eval_hyp2f1(b, _1-c+b, _1-a+b, 1/z)
        #s1 = G(c)*G(b-a)/G(b)/G(c-a) * (-z)**(-a) * h1
        #s2 = G(c)*G(a-b)/G(a)/G(c-b) * (-z)**(-b) * h2
        f1 = gammaprod([c,b-a],[b,c-a])
        f2 = gammaprod([c,a-b],[a,c-b])
        s1 = f1 * (-z)**(_0-a) * h1
        s2 = f2 * (-z)**(_0-b) * h2
        v = s1 + s2
    finally:
        mp.prec = orig
    return +v

#---------------------------------------------------------------------------#
#                      And now the user-friendly versions                   #
#---------------------------------------------------------------------------#

def hyper(as, bs, z):
    """
    Hypergeometric function pFq,

          [ a_1, a_2, ..., a_p |    ]
      pFq [                    |  z ]
          [ b_1, b_2, ..., b_q |    ]

    The parameter lists as and bs may contain real or complex numbers.
    Exact rational parameters can be given as tuples (p, q).
    """
    p = len(as)
    q = len(bs)
    z = convert_lossless(z)
    degree = p, q
    if degree == (0, 1):
        br, bf, bc = parse_param(bs[0])
        if br:
            return sum_hyp0f1_rat(br[0], z)
        return hypsum([], [], [], br, bf, bc, z)
    if degree == (1, 1):
        ar, af, ac = parse_param(as[0])
        br, bf, bc = parse_param(bs[0])
        if ar and br:
            a, b = ar[0], br[0]
            return sum_hyp1f1_rat(a, b, z)
        return hypsum(ar, af, ac, br, bf, bc, z)
    if degree == (2, 1):
        return eval_hyp2f1(as[0],as[1],bs[0],z)
    ars, afs, acs, brs, bfs, bcs = [], [], [], [], [], []
    for a in as:
        r, f, c = parse_param(a)
        ars += r
        afs += f
        acs += c
    for b in bs:
        r, f, c = parse_param(b)
        brs += r
        bfs += f
        bcs += c
    return hypsum(ars, afs, acs, brs, bfs, bcs, z)

def hyp0f1(a, z):
    """Hypergeometric function 0F1. hyp0f1(a,z) is equivalent
    to hyper([], [a], z); see documentation for hyper() for more
    information."""
    return hyper([], [a], z)

def hyp1f1(a,b,z):
    """Hypergeometric function 1F1. hyp1f1(a,b,z) is equivalent
    to hyper([a], [b], z); see documentation for hyper() for more
    information."""
    return hyper([a], [b], z)

def hyp2f1(a,b,c,z):
    """Hypergeometric function 2F1. hyp2f1(a,b,c,z) is equivalent
    to hyper([a,b], [c], z); see documentation for hyper() for more
    information."""
    return hyper([a,b], [c], z)

def funcwrapper(f):
    def g(z):
        orig = mp.prec
        rnd = mp.rounding[0]
        try:
            z = convert_lossless(z)
            mp.prec = orig + 10
            v = f(z)
        finally:
            mp.prec = orig
        return +v
    g.__name__ = f.__name__
    g.__doc__ = f.__doc__
    return g

@extraprec(20, normalize_output=True)
def lower_gamma(a,z):
    """Lower incomplete gamma function gamma(a, z)"""
    z = convert_lossless(z)
    if not isinstance(a, (int, long)):
        a = convert_lossless(a)
    # XXX: may need more precision
    return hyp1f1(1, 1+a, z) * z**a * exp(-z) / a

@extraprec(20, normalize_output=True)
def upper_gamma(a,z):
    """Upper incomplete gamma function Gamma(a, z)"""
    return gamma(a) - lower_gamma(a, z)

@funcwrapper
def erf(z):
    """Error function, erf(z)"""
    return (2/sqrt(pi)*z) * sum_hyp1f1_rat((1,2),(3,2), -z**2)

@funcwrapper
def ellipk(m):
    """Complete elliptic integral of the first kind, K(m). Note that
    the argument is the parameter m = k^2, not the modulus k."""
    if m == 1:
        return inf
    return pi/2 * sum_hyp2f1_rat((1,2),(1,2),(1,1), m)

@funcwrapper
def ellipe(m):
    """Complete elliptic integral of the second kind, E(m). Note that
    the argument is the parameter m = k^2, not the modulus k."""
    if m == 1:
        return m
    return pi/2 * sum_hyp2f1_rat((1,2),(-1,2),(1,1), m)

# TODO: for complex a, b handle the branch cut correctly
@extraprec(15, normalize_output=True)
def agm(a, b):
    """Arithmetic-geometric mean of a and b."""
    a = convert_lossless(a)
    b = convert_lossless(b)
    if not a or not b:
        return a*b
    weps = eps * 16
    half = mpf(0.5)
    while abs(a-b) > weps:
        a, b = (a+b)*half, (a*b)**half
    return a

def jacobi(n, a, b, x):
    """Jacobi polynomial P_n^(a,b)(x)."""
    orig = mp.prec
    try:
        mp.prec = orig + 15
        x = convert_lossless(x)
        v = binomial(n+a,n) * hyp2f1(-n,1+n+a+b,a+1,(1-x)/2)
    finally:
        mp.prec = orig
    return +v

def legendre(n, x):
    """Legendre polynomial P_n(x)."""
    orig = mp.prec
    try:
        mp.prec = orig + 15
        x = convert_lossless(x)
        if not isinstance(n, (int, long)):
            n = convert_lossless(n)
        if x == -1:
            # TODO: hyp2f1 should handle this
            if x == int(x):
                return (-1)**(n + (n>=0)) * mpf(-1)
            return inf
        v = hyp2f1(-n,n+1,1,(1-x)/2)
    finally:
        mp.prec = orig
    return +v

def chebyt(n, x):
    """Chebyshev polynomial of the first kind T_n(x)."""
    orig = mp.prec
    try:
        mp.prec = orig + 15
        x = convert_lossless(x)
        v = hyp2f1(-n,n,0.5,(1-x)/2)
    finally:
        mp.prec = orig
    return +v

def chebyu(n, x):
    """Chebyshev polynomial of the second kind U_n(x)."""
    orig = mp.prec
    try:
        mp.prec = orig + 15
        x = convert_lossless(x)
        v = (n+1) * hyp2f1(-n,n+2,1.5,(1-x)/2)
    finally:
        mp.prec = orig
    return +v

# A Bessel function of the first kind of integer order, J_n(x), is
# given by the power series

#             oo
#             ___         k         2 k + n
#            \        (-1)     / x \
#    J_n(x) = )    ----------- | - |
#            /___  k! (k + n)! \ 2 /
#            k = 0

# Simplifying the quotient between two successive terms gives the
# ratio x^2 / (-4*k*(k+n)). Hence, we only need one full-precision
# multiplication and one division by a small integer per term.
# The complex version is very similar, the only difference being
# that the multiplication is actually 4 multiplies.

# In the general case, we have
# J_v(x) = (x/2)**v / v! * 0F1(v+1, (-1/4)*z**2)

# TODO: for extremely large x, we could use an asymptotic
# trigonometric approximation.

# TODO: recompute at higher precision if the fixed-point mantissa
# is very small

def mpf_jn_series(n, x, prec):
    negate = n < 0 and n & 1
    n = abs(n)
    origprec = prec
    prec += 20 + bitcount(abs(n))
    x = to_fixed(x, prec)
    x2 = (x**2) >> prec
    if not n:
        s = t = MP_ONE << prec
    else:
        s = t = (x**n // int_fac(n)) >> ((n-1)*prec + n)
    k = 1
    while t:
        t = ((t * x2) // (-4*k*(k+n))) >> prec
        s += t
        k += 1
    if negate:
        s = -s
    return make_mpf(from_man_exp(s, -prec, origprec, round_nearest))

def mpc_jn_series(n, z, prec):
    negate = n < 0 and n & 1
    n = abs(n)
    origprec = prec
    prec += 20 + bitcount(abs(n))
    zre, zim = z
    zre = to_fixed(zre, prec)
    zim = to_fixed(zim, prec)
    z2re = (zre**2 - zim**2) >> prec
    z2im = (zre*zim) >> (prec-1)
    if not n:
        sre = tre = MP_ONE << prec
        sim = tim = MP_ZERO
    else:
        re, im = complex_int_pow(zre, zim, n)
        sre = tre = (re // int_fac(n)) >> ((n-1)*prec + n)
        sim = tim = (im // int_fac(n)) >> ((n-1)*prec + n)
    k = 1
    while abs(tre) + abs(tim) > 3:
        p = -4*k*(k+n)
        tre, tim = tre*z2re - tim*z2im, tim*z2re + tre*z2im
        tre = (tre // p) >> prec
        tim = (tim // p) >> prec
        sre += tre
        sim += tim
        k += 1
    if negate:
        sre = -sre
        sim = -sim
    re = from_man_exp(sre, -prec, origprec, round_nearest)
    im = from_man_exp(sim, -prec, origprec, round_nearest)
    return make_mpc((re, im))

def jv(v, x):
    """Bessel function J_v(x)."""
    prec = mp.prec
    x = convert_lossless(x)
    if isinstance(v, int_types):
        if isinstance(x, mpf):
            return mpf_jn_series(v, x._mpf_, prec)
        if isinstance(x, mpc):
            return mpc_jn_series(v, (x.real._mpf_, x.imag._mpf_), prec)
    hx = x/2
    return hx**v * hyp0f1(v+1, -hx**2) / factorial(v)

jn = jv

def j0(x):
    """Bessel function J_0(x)."""
    return jv(0, x)

def j1(x):
    """Bessel function J_1(x)."""
    return jv(1, x)

#---------------------------------------------------------------------------#
#                                                                           #
#                               Miscellaneous                               #
#                                                                           #
#---------------------------------------------------------------------------#


def log_range():
    """Generate log(2), log(3), log(4), ..."""
    prec = mp.prec + 20
    one = 1 << prec
    L = log2_fixed(prec)
    p = 2
    while 1:
        yield mpf((L, -prec))
        s = 0
        u = one
        k = 1
        a = (2*p+1)**2
        while u:
            s += u // k
            u //= a
            k += 2
        L += 2*s//(2*p+1)
        p += 1

@extraprec(30, normalize_output=True)
def lambertw(z, k=0, approx=None):
    """
    lambertw(z,k) gives the kth branch of the Lambert W function W(z),
    defined as the kth solution of z = W(z)*exp(W(z)).

    lambertw(z) == lambertw(z, k=0) gives the principal branch
    value (0th branch solution), which is real for z > -1/e .

    The k = -1 branch is real for -1/e < z < 0. All branches except
    k = 0 have a logarithmic singularity at 0.

    The definition, implementation and choice of branches is based
    on Corless et al, "On the Lambert W function", Adv. Comp. Math. 5
    (1996) 329-359, available online here:
    http://www.apmaths.uwo.ca/~djeffrey/Offprints/W-adv-cm.pdf

    TODO: use a series expansion when extremely close to the branch point
    at -1/e and make sure that the proper branch is chosen there
    """
    z = convert_lossless(z)
    if isnan(z):
        return z
    # We must be extremely careful near the singularities at -1/e and 0
    u = exp(-1)
    if abs(z) <= u:
        if not z:
            # w(0,0) = 0; for all other branches we hit the pole
            if not k:
                return z
            return -inf
        if not k:
            w = z
        # For small real z < 0, the -1 branch behaves roughly like log(-z)
        elif k == -1 and not z.imag and z.real < 0:
            w = log(-z)
        # Use a simple asymptotic approximation.
        else:
            w = log(z)
            # The branches are roughly logarithmic. This approximation
            # gets better for large |k|; need to check that this always
            # works for k ~= -1, 0, 1.
            if k: w += k * 2*pi*j
    else:
        if z == inf: return z
        if z == -inf: return nan
        # Simple asymptotic approximation as above
        w = log(z)
        if k: w += k * 2*pi*j
    # Use Halley iteration to solve w*exp(w) = z
    two = mpf(2)
    weps = ldexp(eps, 15)
    for i in xrange(100):
        ew = exp(w)
        wew = w*ew
        wewz = wew-z
        wn = w - wewz/(wew+ew-(w+two)*wewz/(two*w+two))
        if abs(wn-w) < weps*abs(wn):
            return wn
        else:
            w = wn
    print "Warning: Lambert W iteration failed to converge:", z
    return wn
