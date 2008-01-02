"""
Mathematical constants

"""

import math

from util import *
from floatop import *
from squareroot import *


def constant_memo(f):
    """
    Memoization for calculation of constants using fixed-point
    arithmetic. Only store the value computed with highest precision;
    if a lower or equal precision is requested later, the result can
    be generated directly through downshifting the cached value.
    """
    f.memo_prec = -1
    f.memo_val = None

    def g(prec):
        if prec == f.memo_prec:
            return f.memo_val
        if prec < f.memo_prec:
            return f.memo_val >> (f.memo_prec-prec)
        f.memo_val = f(prec)
        f.memo_prec = prec
        return f.memo_val

    g.__name__ = f.__name__
    g.__doc__ = f.__doc__
    return g

def acot(n, prec, hyperbolic):
    """
    Compute acot of an integer using fixed-point arithmetic. With
    hyperbolic=True, compute acoth. We use the series expansion

                   1        1        1
        acot(n) = ---  -  ----  +  ----  -  ...,
                            3         5
                   n      3 n      5 n

    optimized for integer arguments. In the hyperbolic case, all
    negative terms are changed to positive ones.
    """
    s = t = (1 << prec) // n  # 1 / n
    k = 3
    while 1:
        # Repeatedly divide by k * n**2, and add
        t //= (n*n)
        term = t // k
        if not term:
            break
        # Alternate signs
        if hyperbolic or not k & 2:
            s += term
        else:
            s -= term
        k += 2
    return s

def machin(coefs, prec, hyperbolic=False):
    """
    Evaluate a Machin-like formula, i.e., a linear combination of
    acot(n) or acoth(n) for specific integer values of n, using fixed-
    point arithmetic.

    The input should be a list [(c, n), ...], giving c*acot[h](n) + ...
    """
    extraprec = 10
    s = 0
    for a, b in coefs:
        s += a * acot(b, prec+extraprec, hyperbolic)
    return (s >> extraprec)

def agm_status(prec, step, adiff, verbose_base):
    logdiff = math.log(max(1, adiff), verbose_base)
    digits = int(prec/math.log(verbose_base,2) - logdiff)
    print "  iteration", step, ("(accuracy ~= %i base-%i digits)" % \
       (digits, verbose_base))

def pi_agm(prec, verbose=False, verbose_base=10):
    """
    Compute floor(pi * 2**prec) as a big integer using the Brent-
    Salamin algorithm based on the arithmetic-geometric mean.

    See for example Wikipedia (http://en.wikipedia.org/wiki/Brent-
    Salamin_algorithm) or "Pi and the AGM" by Jonathan and Peter
    Borwein (Wiley, 1987). The algorithm (as stated in the Wikipedia
    article) consists of setting

      a_0 = 1
      b_0 = 1/sqrt(2)
      t_0 = 1/4
      p_0 = 1

    and computing

      a_{n+1} = (a_n + b_n)/2
      b_{n+1} = sqrt(a_n * b_n)
      t_{n+1} = t_n - p_n*(a_n - a_{n+1})**2
      p_{n+1} = 2*p_n

    for n = 0, 1, 2, 3, ..., after which the approximation is given by
    pi ~= (a_n + b_n)**2 / (4*t_n). Each step roughly doubles the
    number of correct digits.
    """
    extraprec = 50
    prec += extraprec

    # Initialial values. a, b and t are fixed-point numbers
    a = 1 << prec
    b = sqrt_fixed2(a >> 1, prec)
    t = a >> 2
    p = 1

    step = 1
    while 1:
        an = (a + b) >> 1
        adiff = a - an
        if verbose:
            agm_status(prec, step, adiff, verbose_base)
        # No change in a
        if p > 16 and abs(adiff) < 1000:
            break
        prod = (a * b) >> prec
        b = sqrt_fixed2(prod, prec)
        t = t - p*((adiff**2) >> prec)
        p = 2*p
        a = an
        step += 1
    if verbose:
        print "  final division"
    pi = ((((a+b)**2) >> 2) // t)
    return pi >> extraprec

@constant_memo
def pi_fixed(prec):
    """
    Compute floor(pi * 2**prec) as a big integer.

    For low precisions, Machin's formula pi = 16*acot(5)-4*acot(239)
    is used. For high precisions, the more efficient arithmetic-
    geometric mean iteration is used.
    """
    if prec < 10000:
        return machin([(16, 5), (-4, 239)], prec)
    else:
        return pi_agm(prec)

def fpi(prec, rounding):
    """Compute a floating-point approximation of pi"""
    return normalize(pi_fixed(prec+5), -prec-5, prec, rounding)


# Logarithms of integers can be computed easily using
# Machin-like formulas

@constant_memo
def log2_fixed(prec):
    return machin([(18, 26), (-2, 4801), (8, 8749)], prec, True)

def flog2(prec, rounding):
    return normalize(log2_fixed(prec+5), -prec-5, prec, rounding)


@constant_memo
def log10_fixed(prec):
    return machin([(46, 31), (34, 49), (20, 161)], prec, True)

def flog10(prec, rounding):
    return normalize(log10_fixed(prec+5), -prec-5, prec, rounding)


"""
Euler's constant (gamma) is computed using the Brent-McMillan formula,
gamma ~= A(n)/B(n) - log(n), where

  A(n) = sum_{k=0,1,2,...} (n**k / k!)**2 * H(k)
  B(n) = sum_{k=0,1,2,...} (n**k / k!)**2
  H(k) = 1 + 1/2 + 1/3 + ... + 1/k

The error is bounded by O(exp(-4n)). Choosing n to be a power
of two, 2**p, the logarithm becomes particularly easy to calculate.

Reference:
Xavier Gourdon & Pascal Sebah, The Euler constant: gamma
http://numbers.computation.free.fr/Constants/Gamma/gamma.pdf
"""

@constant_memo
def gamma_fixed(prec):
    prec += 30
    # choose p such that exp(-4*(2**p)) < 2**-n
    p = int(math.log((prec/4) * math.log(2), 2)) + 1
    n = 1<<p
    r = one = 1<<prec
    H, A, B, npow, k, d = 0, 0, 0, 1, 1, 1
    while r:
        A += (r * H) >> prec
        B += r
        r = r * (n*n) // (k*k)
        H += one // k
        k += 1
    S = ((A<<prec) // B) - p*log2_fixed(prec)
    return S >> 30

def fgamma(prec, rounding):
    return normalize(gamma_fixed(prec+5), -prec-5, prec, rounding)
