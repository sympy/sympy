"""
Commonly used mathematical constants (pi, log(2), log(10), gamma, ...)

There are two functions for each constant: const_fixed(prec) returns
a big integer n such that const ~= n * 2**-prec, and const_float()
returns Float(const), rounded to the current Float working precision
"""

from math import log as _clog
from float_ import Float
from utils_ import global_options

# Only re-compute a constant if the precision level is raised
def _constmemo(f):
    f.memo_prec = -1
    f.memo_val = None
    def calc(prec):
        if prec == f.memo_prec: return f.memo_val
        if prec < f.memo_prec: return f.memo_val >> (f.memo_prec-prec)
        f.memo_val = f(prec)
        f.memo_prec = prec
        return f.memo_val
    return calc

# Evaluate a Machin-like formula, i.e., a rational combination of
# of acot(n) or acoth(n) for specific integer values of n
def _machin(coefs, prec, hyperbolic=False):
    prec += 10
    def acot(x):
        # Use the basic series expansion for atan/acot, optimized
        # for integer arguments
        s = w = (1<<prec)//x; x2 = x*x;  n = 3
        while 1:
            w //= x2
            term = w // n
            if not term: break
            if hyperbolic or n & 2 == 0: s += term
            else: s -= term
            n += 2
        return s
    s = 0
    for a, b in coefs:
        s += a * acot(b)
    return s >> 10

#----------------------------------------------------------------------
# Pi
#

"""
At low precision, pi can be calculated easily using Machin's formula
pi = 16*acot(5)-4*acot(239).

For high precision, we use the Brent-Salamin algorithm based on the
arithmetic-geometric mean. See for example Wikipedia
(http://en.wikipedia.org/wiki/Brent-Salamin_algorithm)
or "Pi and the AGM" by Jonathan and Peter Borwein (Wiley, 1987).
The algorithm (as stated in the Wikipedia article) consists of
setting

  a_0 = 1;  b_0 = 1/sqrt(2);  t_0 = 1/4;  p_0 = 1

and computing

  a_{n+1} = (a_n + b_n)/2
  b_{n+1} = sqrt(a_n * b_n)
  t_{n+1} = t_n - p_n*(a_n - a_{n+1})**2
  p_{n+1} = 2*p_n

for n = 0, 1, 2, 3, ..., after which the approximation is given by

  pi ~= (a_n + b_n)**2 / (4*t_n).

Each step roughly doubles the number of correct digits.
"""

def _pi_agm(prec):
    from functions import _sqrt_fixed2
    prec += 50
    a = 1 << prec
    if "verbose" in global_options:
        print "  computing initial square root..."
    b = _sqrt_fixed2(a >> 1, prec)
    t = a >> 2
    p = 1
    step = 1
    while 1:
        an = (a+b)>>1
        adiff = a - an
        if "verbose" in global_options:
            base = global_options.get("verbose_base", 10)
            try:
                logdiff = _clog(adiff, base)
            except ValueError:
                logdiff = 0
            digits = int(prec/_clog(base,2) - logdiff)
            print "  iteration", step, ("(accuracy ~= %i base-%i digits)" % (digits, base))
        if p > 16 and abs(adiff) < 1000:
            break
        prod = (a*b)>>prec
        b = _sqrt_fixed2(prod, prec)
        t = t - p*((adiff**2) >> prec)
        p = 2*p
        a = an
        step += 1
    if "verbose" in global_options:
        print "  final division"
    return ((((a+b)**2) >> 2) // t) >> 50

@_constmemo
def pi_fixed(prec):
    if prec < 10000:
        return _machin([(16, 5), (-4, 239)], prec)
    else:
        return _pi_agm(prec)

def pi_float():
    prec = Float._prec + 5
    return Float((pi_fixed(prec), -prec))

#----------------------------------------------------------------------
# Logarithms
#

# Logarithms of integers can be computed easily using
# Machin-like formulas

@_constmemo
def log2_fixed(prec):
    return _machin([(18, 26), (-2, 4801), (8, 8749)], prec, True)

def log2_float():
    prec = Float._prec + 5
    return Float((log2_fixed(prec), -prec))

@_constmemo
def log10_fixed(prec):
    return _machin([(46, 31), (34, 49), (20, 161)], prec, True)

def log10_float():
    prec = Float._prec + 5
    return Float((log10_fixed(prec), -prec))

#----------------------------------------------------------------------
# Other constants
#

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

@_constmemo
def gamma_fixed(prec):
    # XXX: may need even more extra precision
    prec += 30
    # choose p such that exp(-4*(2**p)) < 2**-n
    p = int(_clog((prec/4) * _clog(2), 2)) + 1
    n = 1<<p; r=one=1<<prec
    H, A, B, npow, k, d = 0, 0, 0, 1, 1, 1
    while r:
        A += (r * H) >> prec
        B += r
        r = r * (n*n) // (k*k)
        H += one // k
        k += 1
    S = ((A<<prec) // B) - p*log2_fixed(prec)
    return S >> 30

def gamma_float():
    prec = Float._prec + 5
    return Float((gamma_fixed(prec), -prec))
