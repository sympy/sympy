"""
Numerical implementations of special functions (gamma, ...)
"""

from float_ import Float
from constants import pi_float
from functions import exp, log, sqrt, sin
from utils_ import make_fixed


#----------------------------------------------------------------------
# Gamma function
#

"""
We compute the gamma function using Spouge's approximation

    x! = (x+a)**(x+1/2) * exp(-x-a) * [c_0 + S(x) + eps]

where S(x) is the sum c_k/(x+k) from k = 1 to a-1 and the coefficients
are given by

  c_0 = sqrt(2*pi)

         (-1)**(k-1)
  c_k =  ----------- (a-k)**(k-1/2) exp(-k+a),  k = 1,2,...,a-1
          (k - 1)!

If we choose a = int(1.26*n), the error eps is less than 10**-n
for any x in the right complex half-plane (assuming a > 2). In practice,
it seems that a can be chosen quite a bit lower (30-50%); this needs
to be looked into.

Reference:
John L. Spouge, "Computation of the gamma, digamma, and trigamma
functions", SIAM Journal on Numerical Analysis 31 (1994), no. 3, 931-944.
"""

_spouge_cache = {}

# Calculate coefficients for Spouge's approximation
def _calc_spouge(a, prec):
    c = [0] * a
    # We'll store the coefficients as fixed-point numbers but calculate
    # them as Floats for convenience. The initial terms are huge, so we
    # need to allocate extra bits to ensure full accuracy. The integer
    # part of the largest term has ~= exp(a) or 2**(1.4*a) bits.
    floatprec = prec + int(a*1.4)
    Float.store()
    Float.setprec(floatprec)
    b = exp(a-1)
    e = exp(1)
    c[0] = make_fixed(sqrt(2*pi_float()), prec)
    for k in range(1, a):
        t = ((-1)**(k-1) * (a-k)**k)
        u = t * b / sqrt(a-k)
        c[k] = make_fixed(u, prec)
        b = b / (e * k)
        # print "computing coefficients: %02f" % (100.0 * k / a), "% done"
    Float.revert()
    return c

# Cached lookup of coefficients
def _get_spouge(prec):
    if prec in _spouge_cache:
        return _spouge_cache[prec]
    for p in _spouge_cache:
        # Coefficients calculated for a slightly higher precision are ok too
        if 0.8 <= float(p)/prec < 1:
            return _spouge_cache[p]
    a = max(3, int(0.35*prec))  # ~= 1.26*n
    coefs = _calc_spouge(a, prec)
    _spouge_cache[prec] = (prec, a, coefs)
    return _spouge_cache[prec]

def gamma(x):
    """
    gamma(x) -- calculate the gamma function of the real number x

    x must not be a negative integer or 0
    """
    if not isinstance(x, Float):
        x = Float(x)
    if x < 0.25:
        if x == int(x):
            raise ZeroDivisionError, "gamma function pole"
        # Use reflection formula
        Float._prec += 3
        p = pi_float()
        g = p / (sin(p*x) * gamma(1-x))
        Float._prec -= 3
        return +g
    x -= 1
    prec, a, c = _get_spouge(Float.getprec() + 10)
    # Perform summation in fixed-point mode
    xf = make_fixed(x, prec)
    s = c[0]
    for k in xrange(1, a):
        term = (c[k] << prec) // (xf + (k << prec))
        s += term
    Float._prec += 8
    g = exp(log(x+a)*(x+Float(0.5))) * exp(-(x+a)) * Float((s, -prec))
    Float._prec -= 8
    return +g
