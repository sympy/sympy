from __future__ import print_function, division

from sympy.mpmath.libmp import (fzero,
    from_man_exp, from_int, from_rational,
    fone, fhalf, bitcount, to_int, to_str, mpf_mul, mpf_div, mpf_sub,
    mpf_add, mpf_sqrt, mpf_pi, mpf_cosh_sinh, pi_fixed, mpf_cos)
from sympy.core.numbers import igcd
import math
from sympy.core.compatibility import xrange


def _a(n, j, prec):
    """Compute the inner sum in the HRR formula."""
    if j == 1:
        return fone
    s = fzero
    pi = pi_fixed(prec)
    for h in xrange(1, j):
        if igcd(h, j) != 1:
            continue
        # & with mask to compute fractional part of fixed-point number
        one = 1 << prec
        onemask = one - 1
        half = one >> 1
        g = 0
        if j >= 3:
            for k in xrange(1, j):
                t = h*k*one//j
                if t > 0:
                    frac = t & onemask
                else:
                    frac = -((-t) & onemask)
                g += k*(frac - half)
        g = ((g - 2*h*n*one)*pi//j) >> prec
        s = mpf_add(s, mpf_cos(from_man_exp(g, -prec), prec), prec)
    return s


def _d(n, j, prec, sq23pi, sqrt8):
    """
    Compute the sinh term in the outer sum of the HRR formula.
    The constants sqrt(2/3*pi) and sqrt(8) must be precomputed.
    """
    j = from_int(j)
    pi = mpf_pi(prec)
    a = mpf_div(sq23pi, j, prec)
    b = mpf_sub(from_int(n), from_rational(1, 24, prec), prec)
    c = mpf_sqrt(b, prec)
    ch, sh = mpf_cosh_sinh(mpf_mul(a, c), prec)
    D = mpf_div(mpf_sqrt(j, prec), mpf_mul(mpf_mul(sqrt8, b), pi), prec)
    E = mpf_sub(mpf_mul(a, ch), mpf_div(sh, c, prec), prec)
    return mpf_mul(D, E)


def npartitions(n, verbose=False):
    """
    Calculate the partition function P(n), i.e. the number of ways that
    n can be written as a sum of positive integers.

    P(n) is computed using the Hardy-Ramanujan-Rademacher formula,
    described e.g. at http://mathworld.wolfram.com/PartitionFunctionP.html

    The correctness of this implementation has been tested for 10**n
    up to n = 8.

    Examples
    ========

    >>> from sympy.ntheory import npartitions
    >>> npartitions(25)
    1958
    """
    n = int(n)
    if n < 0:
        return 0
    if n <= 5:
        return [1, 1, 2, 3, 5, 7][n]
    # Estimate number of bits in p(n). This formula could be tidied
    pbits = int((math.pi*(2*n/3.)**0.5 - math.log(4*n))/math.log(10) + 1) * \
        math.log(10, 2)
    prec = p = int(pbits*1.1 + 100)
    s = fzero
    M = max(6, int(0.24*n**0.5 + 4))
    sq23pi = mpf_mul(mpf_sqrt(from_rational(2, 3, p), p), mpf_pi(p), p)
    sqrt8 = mpf_sqrt(from_int(8), p)
    for q in xrange(1, M):
        a = _a(n, q, p)
        d = _d(n, q, p, sq23pi, sqrt8)
        s = mpf_add(s, mpf_mul(a, d), prec)
        if verbose:
            print("step", q, "of", M, to_str(a, 10), to_str(d, 10))
        # On average, the terms decrease rapidly in magnitude. Dynamically
        # reducing the precision greatly improves performance.
        p = bitcount(abs(to_int(d))) + 50
    return int(to_int(mpf_add(s, fhalf, prec)))

__all__ = ['npartitions']
