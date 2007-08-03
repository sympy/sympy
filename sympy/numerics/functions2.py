"""
Numerical implementations of special functions (gamma, ...)
"""

from float_ import Float, ComplexFloat
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


# This function computes S
def _spouge_sum(x, prec, a, c):
    if isinstance(x, Float):
        # Regular fixed-point summation
        x = make_fixed(x, prec)
        s = c[0]
        for k in xrange(1, a):
            s += (c[k] << prec) // (x + (k << prec))
        return Float((s, -prec))
    elif isinstance(x, (Rational, int, long)):
        # Here we can save some work
        if isinstance(x, (int, long)):
            p, q = x, 1
        else:
            p, q = x.p, x.q
        s = c[0]
        for k in xrange(1, a):
            s += c[k] * q // (p+q*k)
        return Float((s, -prec))
    elif isinstance(x, ComplexFloat):
        """
        For a complex number a + b*I, we have

              c_k          (a+k)*c_k     b * c_k
        -------------  =   ---------  -  ------- * I
        (a + b*I) + k          M            M

                       2    2      2   2              2
        where M = (a+k)  + b   = (a + b ) + (2*a*k + k )
        """
        re = make_fixed(x.real, prec)
        im = make_fixed(x.imag, prec)
        sre, sim = c[0], 0
        mag = ((re**2)>>prec) + ((im**2)>>prec)
        for k in xrange(1, a):
            M = mag + re*(2*k) + ((k**2) << prec)
            sre += (c[k] * (re + (k << prec))) // M
            sim -= (c[k] * im) // M
        return ComplexFloat(Float((sre, -prec)), Float((sim, -prec)))


def gamma(x):
    """
    gamma(x) -- calculate the gamma function of a real or complex
    number x.
    
    x must not be a negative integer or 0
    """
    Float.store()
    Float._prec += 2

    if isinstance(x, complex):
        x = ComplexFloat(x)
    elif not isinstance(x, (Float, ComplexFloat, Rational, int, long)):
        x = Float(x)

    if isinstance(x, (ComplexFloat, complex)):
        re, im = x.real, x.imag
    else:
        re, im = x, 0

    # For negative x (or positive x close to the pole at x = 0),
    # we use the reflection formula
    if re < 0.25:
        if re == int(re) and im == 0:
            raise ZeroDivisionError, "gamma function pole"
        Float._prec += 3
        p = pi_float()
        g = p / (sin(p*x) * gamma(1-x))
    else:
        x -= 1
        prec, a, c = _get_spouge_coefficients(Float.getprec()+7)
        s = _spouge_sum(x, prec, a, c)
        if not isinstance(x, (Float, ComplexFloat)):
            x = Float(x)
        # TODO: higher precision may be needed here when the precision
        # and/or size of x are extremely large
        Float._prec += 10
        g = exp(log(x+a)*(x+Float(0.5))) * exp(-x-a) * s

    Float.revert()
    return +g
