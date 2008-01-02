"""
This module provides functions for computing approximate square roots
(in fixed-point or floating-point form) of positive real numbers.
"""

from util import *
from floatop import *


def sqrt_initial(y, prec):
    """Given y = floor(x * 2**prec), compute floor(sqrt(x) * 2**50),
    i.e. calculate a 50-bit approximation of the square root. This is
    done quickly using regular floating-point arithmetic. It is
    assumed that x ~= 1."""

    # Two cases; the second avoids overflow
    if prec < 200: return int(y**0.5 * 2.0**(50 - prec*0.5))
    else:          return int((y >> (prec-100))**0.5)


# XXX: doesn't work
def invsqrt_initial(y, prec):
    """Like sqrt_initial, but computes 1/sqrt(y) instead of sqrt(y)."""
    if prec < 200: return int(y**-0.5 * 2.0**(50 + prec*0.5))
    else:          return int((y >> (prec-100)) ** -0.5)


def sqrt_fixed(y, prec):
    """
    Square root of a fixed-point number.

    Given the big integer y = floor(x * 2**prec), this function returns
    floor(r * 2**prec) where r = sqrt(x).

    We start from a 50-bit estimate for r generated with ordinary
    floating-point arithmetic, and then refines the value to full
    accuracy using the iteration

                 1  /        y  \
        r     = --- | r  +  --- |
         n+1     2  \ n     r_n /

    which is simply Newton's method applied to the equation r**2 = y.

    Newton's method doubles the accuracy with each step. We make use of
    this fact by only using a precision corresponding to the current
    accuracy during intermediate iterations. For example, with a 50-bit
    accurate r_1, r_2 can be computed using 100-bit precision, r_3
    using 200-bit precision, and so on. (In practice, the precision
    levels must be chosen slightly more conservatively to account for
    rounding errors in the last one or two bits.)

    It is assumed that x ~= 1 (the main fsqrt() function fiddles with
    the exponent of the input to reduce it to unit magnitude before
    passing it here.)

    TODO: it would be possible to specify separate precision levels
    for the input and output, which could be useful when calculating
    pure-integer square roots.
    """

    r = sqrt_initial(y, prec)
    extra = 10
    prevp = 50

    for p in giant_steps(50, prec+extra):

        # Explanation: in the first term, we shift by the appropriate number
        # of bits to convert r from the previous precision to the current one.
        # The "-1" divides by two in the same step.

        # In the second term, we do a fixed-point division using the usual
        # formula (y<<precision)//r. The precision term is a bit messy and
        # takes into account the fact that y, r_n and r_{n+1} all have
        # different precision levels. As before, the "-1" divides by two.
        r = lshift_quick(r, p-prevp-1) + (lshift_quick(y, p+prevp-prec-1)//r)

        prevp = p

    return r >> extra


def sqrt_fixed2(y, prec):
    """
    This function is essentially equivalent to sqrt_fixed (see its
    documentation), but uses an asymptotically faster algorithm.

    Instead of using Newton's method to calculate sqrt(y) directly,
    we calculate 1/sqrt(y) with Newton's method and multiply by y to
    obtain sqrt(y). The Newton iteration for 1/sqrt(y) is

                 r
                  n      /            2 \
        r    =  ----  *  | 3  - y * r   |.
         n+1      2      \           n  /

    This is slightly slower at low precision levels since it requires
    three multiplications in each step, as opposed to the single
    division in the Newton iteration for sqrt(y).

    However, since Python uses Karatsuba algorithm for multiplication,
    three multiplications can be performed much more quickly than a
    single division at high precision levels. In practice, the cutoff
    where sqrt_fixed2 becomes faster than sqrt_fixed seems to be around
    60,000 bits.
    """

    # XXX
    from mpmath.lib import to_float, normalize, ROUND_FLOOR
    r = to_float(normalize(y, -prec, 64, ROUND_FLOOR)) ** -0.5
    r = int(r * 2**50)

    # r = invsqrt_initial(y, prec)

    extra = 10
    prevp = 50

    for p in giant_steps(50, prec+extra):

        # This is even messier than in sqrt_fixed. As many shifts as possible
        # have been combined together for optimal speed, at a slight expense
        # of legibility.

        # Compute r**2 at precision p.
        r2 = rshift_quick(r*r, 2*prevp - p)

        # A = r, converted from precision prevp to p
        A = lshift_quick(r, p-prevp)

        # S = y * r2, computed at precision p. We shift y by '-prec' to
        # account for its initial precision, and by 'p' for the fixed-point
        # multiplication
        S = (lshift_quick(y, p-prec) * r2) >> p

        # B = (3-S) and finally the outer product, both done at precision p
        B = (3<<p) - S
        r = (A*B) >> (p+1)

        prevp = p

    # Finally, sqrt(y) = y * (1/sqrt(y))
    r = (r * y) >> prec

    return r >> extra




def fsqrt(s, prec, rounding):
    """
    Floating-point square root.

    Returns a tuple representing the square root of s, rounded to the
    nearest floating-point number in the specified rounding direction.
    The input must be a tuple representing a nonnegative floating-point
    number.
    """
    if s == fone:
        return fone

    man, exp, bc = s

    if not man:
        return fzero

    # Convert to a fixed-point number with prec2 bits. Adjust
    # exponents to be even so that they can be divided in half
    prec2 = prec + 10 + (prec & 1)

    if exp & 1:
        exp -= 1
        man <<= 1
        bc += 1

    # Mantissa may have more bits than we need. Trim it down.
    shift = bc - prec2
    shift -= shift & 1
    man = rshift_quick(man, shift)

    if prec < 65000:
        man = sqrt_fixed(man, prec2)
    else:
        man = sqrt_fixed2(man, prec2)

    return normalize(man, (exp+shift-prec2)>>1, prec, rounding)



def fhypot(x, y, prec, rounding):
    if y == fzero:
        return fabs(x, prec, rounding)
    if x == fzero:
        return fabs(y, prec, rounding)

    RF = ROUND_FLOOR
    hypot2 = fadd(fmul(x,x,prec+4,RF), fmul(y,y,prec+4,RF), prec+4, RF)

    return fsqrt(hypot2, prec, rounding)

