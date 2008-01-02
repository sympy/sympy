"""
Arithmetic operations and functions on complex numbers, represented on
rectangular form as tuples of real floating-point numbers.

This module is quite compact, since most of the dirty work can be
delegated to functions that work on real numbers.
"""

from util import *
from floatop import *
from functions import *
from squareroot import *


# Use fastest rounding mode for intermediate calculations
RF = ROUND_FLOOR


#----------------------------------------------------------------------
# Complex parts
#

def fcabs(a, b, prec, rounding):
    """
    Absolute value of a complex number, |a+bi|. Returns a single
    real number.
    """
    return fhypot(a, b, prec, rounding)


#----------------------------------------------------------------------
# Arithmetic
#

def fcmul(a, b, c, d, prec, rounding):
    """
    Complex multiplication.

    Returns the real and imaginary part of (a+bi)*(c+di), rounded to
    the specified precision. The rounding mode applies to the real and
    imaginary parts separately.

    Implemented the straightforward way. TODO: use a more stable
    algorithm to avoid cancellation.
    """

    # All-real case
    if b == d == fzero:
        return fmul(a, c, prec, rounding), fzero

    ep = prec + 10

    re = fsub(fmul(a,c, ep, RF), fmul(b,d, ep, RF), prec, rounding)
    im = fadd(fmul(a,d, ep, RF), fmul(b,c, ep, RF), prec, rounding)
    return re, im


#----------------------------------------------------------------------
# Square root
#

def fcsqrt(a, b, prec, rounding):
    """
    Complex square root (principal branch).

    We have sqrt(a+bi) = sqrt((r+a)/2) + b/sqrt(2*(r+a))*i where
    r = abs(a+bi), when a+bi is not a negative real number.
    """

    if a == b == fzero:
        return (a, b)

    # When a+bi is a negative real number, we get a real sqrt times i
    if a[0] < 0 and b == fzero:
        im = fsqrt(fneg_exact(a), prec, rounding)
        return (fzero, im)

    ep = prec+20

    t  = fadd(fcabs(a, b, ep, RF), a, ep, RF)  # t = abs(a+bi) + a
    u  = fmul(t, fhalf, ep, RF)                # u = t / 2
    re = fsqrt(u, prec, rounding)              # re = sqrt(u)
    v  = fmul(t, ftwo, ep, RF)                 # v = t * 2
    w  = fsqrt(v, ep, RF)                      # w = sqrt(v)
    im = fdiv(b, w, prec, rounding)            # im = b / w

    return re, im


#----------------------------------------------------------------------
# Exp and log
#

def fcexp(a, b, prec, rounding):
    """
    Complex exponential function.

    We use the direct formula exp(a+bi) = exp(a) * (cos(b) + sin(b)*i)
    for the computation. This formula is very nice because it is
    perfectly stable; since we just do real multiplications, the only
    numerical errors that can creep in are single-ulp rounding errors.

    The formula is efficient since mpmath's real exp is quite fast and
    since we can compute cos and sin simultaneously.

    It is no problem if a and b are large; if the implementations of
    exp/cos/sin are accurate and efficient for all real numbers, then
    so is this function for all complex numbers.
    """

    mag = fexp(a, prec+4, rounding)

    c, s = cos_sin(b, prec+4, rounding)

    re = fmul(mag, c, prec, rounding)
    im = fmul(mag, s, prec, rounding)
    return re, im

# TODO: log


#----------------------------------------------------------------------
# Trigonometric functions
#

def fccos(a, b, prec, rounding):
    """
    Complex cosine.

    The formula used is cos(a+bi) = cos(a)*cosh(b) - sin(a)*sinh(b)*i.

    The same comments apply as for the complex exp: only real
    multiplications are performed, so no cancellation errors are
    possible. The formula is also efficient since we can compute both
    pairs (cos, sin) and (cosh, sinh) in single steps.
    """
    ep = prec + 6

    c, s = cos_sin(a, ep, RF)
    ch, sh = cosh_sinh(b, ep, RF)

    re = fmul(c, ch, prec, rounding)
    im = fmul(s, sh, prec, rounding)

    return re, fneg_exact(im)


def fcsin(a, b, prec, rounding):
    """
    Complex sine.

    We have sin(a+bi) = sin(a)*cosh(b) + cos(a)*sinh(b)*i.
    See the docstring for fccos for additional comments.
    """
    ep = prec + 6

    c, s = cos_sin(a, ep, RF)
    ch, sh = cosh_sinh(b, ep, RF)

    re = fmul(s, ch, prec, rounding)
    im = fmul(c, sh, prec, rounding)

    return re, im


# TODO: complex tan


#----------------------------------------------------------------------
# Hyperbolic functions
#

def fccosh(a, b, prec, rounding):
    """Complex hyperbolic cosine. Computed as cosh(z) = cos(z*i)."""
    return fccos(b, fneg_exact(a), prec, rounding)


def fcsinh(a, b, prec, rounding):
    """Complex hyperbolic sine. Computed as sinh(z) = -i*sin(z*i)."""
    b, a = fcsin(b, a, prec, rounding)
    return a, b
