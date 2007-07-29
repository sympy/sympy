"""
Numerical implementations of special functions (gamma, ...)
"""

from float_ import Float
from constants import pi_float
from functions import exp, log, sqrt, sin
from utils_ import make_fixed

from sympy import Rational


#---------------------------------------------------------------------------#
#                                                                           #
#                              Gamma function                               #
#                                                                           #
#---------------------------------------------------------------------------#

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


#----------------------------------------------------------------------
#
# We first implement a helper function for calculating the coefficients
# c_k and caching them so that they can be re-used for multiple gamma
# function evaluations
#

_spouge_cache = {}

def _calc_spouge_coefficients(a, prec):
    """
    Calculate Spouge coefficients for approximation with parameter a.
    Return a list of big integers representing the coefficients in
    fixed-point form with a precision of prec bits.
    """

    # We'll store the coefficients as fixed-point numbers but calculate
    # them as Floats for convenience. The initial terms are huge, so we
    # need to allocate extra bits to ensure full accuracy. The integer
    # part of the largest term has size ~= exp(a) or 2**(1.4*a)
    floatprec = prec + int(a*1.4)
    Float.store()
    Float.setprec(floatprec)

    c = [0] * a
    b = exp(a-1)
    e = exp(1)
    c[0] = make_fixed(sqrt(2*pi_float()), prec)
    for k in range(1, a):
        # print "%02f" % (100.0 * k / a), "% done"
        c[k] = make_fixed(((-1)**(k-1) * (a-k)**k) * b / sqrt(a-k), prec)
        # Divide off e and k instead of computing exp and k! from scratch
        b = b / (e * k)

    Float.revert()
    return c

# Cached lookup of coefficients
def _get_spouge_coefficients(prec):

    # This exact precision has been used before
    if prec in _spouge_cache:
        return _spouge_cache[prec]

    for p in _spouge_cache:
        # Coefficients calculated for a slightly higher precision are ok
        # too. But if the difference is too big, we're better off
        # starting from scratch
        if 0.8 <= float(p)/prec < 1:
            return _spouge_cache[p]

    # Here we estimate the value of a based on Spouge's inequality for
    # the relative error
    a = max(3, int(0.39*prec))  # ~= 1.26*n

    # Compute and return
    coefs = _calc_spouge_coefficients(a, prec)
    _spouge_cache[prec] = (prec, a, coefs)
    return _spouge_cache[prec]


#----------------------------------------------------------------------
# Helper functions for computing the sum part, S(x)
#

# Optimizing for different types of numbers makes the computation
# several times faster in some common cases

def _spouge_sum_float(x, prec, a, c):
    xf = make_fixed(x, prec)
    s = c[0]
    for k in xrange(1, a):
        s += (c[k] << prec) // (xf + (k << prec))
    return Float((s, -prec))

def _spouge_sum_rational(p, q, prec, a, c):
    s = c[0]
    for k in xrange(1, a):
        s += c[k] * q // (p+q*k)
    return Float((s, -prec))

# to be implemented
def _spouge_sum_complex():
    pass


#----------------------------------------------------------------------
# Main function
#

def gamma(x):
    """
    gamma(x) -- calculate the gamma function of the real number x.

    x must not be a negative integer or 0
    """
    Float.store()
    Float._prec += 2

    if not isinstance(x, (Float, Rational, int, long)):
        x = Float(x)

    # For negative x (or positive x close to the pole at x = 0),
    # we use the reflection formula
    if x < 0.25:
        if x == int(x):
            raise ZeroDivisionError, "gamma function pole"
        Float._prec += 3
        p = pi_float()
        g = p / (sin(p*x) * gamma(1-x))
        Float._prec -= 3
        return +g

    x -= 1
    prec, a, c = _get_spouge_coefficients(Float.getprec()+7)
    # Sum the series
    if isinstance(x, (int, long)):
        x = Rational(x)
    if isinstance(x, Rational):
        s = _spouge_sum_rational(x.p, x.q, prec, a, c)
        x = Float(x)
    else:
        x = Float(x)
        s = _spouge_sum_float(x, prec, a, c)

    # TODO: higher precision may be needed here when the precision
    # and/or size of x are extremely large
    Float._prec += 10
    g = exp(log(x+a)*(x+Float(0.5))) * exp(-x-a) * s

    Float.revert()
    return +g
