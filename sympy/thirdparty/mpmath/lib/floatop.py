"""
Functions for basic operations on raw mpfs: normalization, comparison,
addition, subtraction, multiplication, division, integer powers.
"""

from util import *
import random as _random

# Some commonly needed raw mpfs
fzero = (0, 0, 0)
fone = (1, 0, 1)
ftwo = (1, 1, 1)
ften = (5, 1, 3)
fhalf = (1, -1, 1)


#-----------------------------------------------------------------------------
#
# normalize() is the workhorse function in mpmath. All floating-point
# operations are implemented according to the pattern
#
#  1) convert floating-point problem to an equivalent integer problem
#  2) solve integer problem using Python integer arithmetic
#  3) use normalize() to convert the integer solution to a canonical
#     floating-point number
#
# A number of hacks are used in normalize() to reduce the overhead of step
# (3) as far as possible.
#

# Pre-computing and avoiding calls to trailing_zeros() in
# normalize improves performance at 15-digit precision by ~15%
shift_table = map(trailing_zeros, range(256))

# The general method for counting bits is to use the bitcount() function,
# which calls math.log. However, for prec ~< 200 bits, converting to a 
# hex string and using table lookup for the last four bits is slightly
# faster. This idea was taken from decimal.py
correction = {
        '0': 4, '1': 3, '2': 2, '3': 2,
        '4': 1, '5': 1, '6': 1, '7': 1,
        '8': 0, '9': 0, 'a': 0, 'b': 0,
        'c': 0, 'd': 0, 'e': 0, 'f': 0}


def normalize(man, exp, prec, rounding):
    """
    normalize(man, exp, prec, rounding) -> return tuple representing
    a fully rounded and normalized raw mpf with value (man * 2**exp)

    The mantissa is rounded in the specified direction if the number of
    exceeds the precision. Trailing zero bits are also stripped from
    the mantissa to ensure that the representation is canonical.
    """

    # The bit-level operations below assume a nonzero mantissa
    if not man:
        return fzero
    # Count bits
    if prec < 200:
        # Inlined for a small speed boost at low precision
        hex_n = "%x" % abs(man)
        bc = 4*len(hex_n) - correction[hex_n[0]]
    else:
        bc = bitcount(man)
    # Cut mantissa down to size
    if bc > prec:
        man = rshift(man, bc-prec, rounding)
        exp += (bc - prec)
        bc = prec
    # Strip trailing zeros
    if not man & 1:
        # To within the nearest byte
        while not man & 0xff:
            man >>= 8
            exp += 8
            bc -= 8
        t = shift_table[man & 0xff]
        man >>= t
        exp += t
        bc -= t
    # If result is +/- a power of two due to rounding up in rshift(),
    # bc may be wrong
    if man == 1 or man == -1:
        bc = 1
    return (man, exp, bc)


# Equivalent to normalize(), but takes a full (man, exp, bc) tuple
def fpos(s, prec, rounding):
    """Calculate 0+s for a raw mpf (i.e., just round s to the specified
    precision, or return s unchanged if its mantissa is smaller than
    the precision)."""
    return normalize(s[0], s[1], prec, rounding)


#-----------------------------------------------------------------------------
# Comparison operations
#

def feq(s, t):
    """Test equality of two raw mpfs. (This is simply tuple comparion;
    this function is provided only for completeness)."""
    return s == t


def fcmp(s, t):
    """Compare the raw mpfs s and t. Return -1 if s < t, 0 if s == t,
    and 1 if s > t. (Same convention as Python's cmp() function.)"""

    # In principle, a comparison amounts to determining the sign of s-t.
    # A full subtraction is relatively slow, however, so we first try to
    # look at the components.
    sman, sexp, sbc = s
    tman, texp, tbc = t

    # Very easy cases: check for zeros and opposite signs
    if not tman: return cmp(sman, 0)
    if not sman: return cmp(0, tman)
    if sman > 0 and tman < 0: return 1
    if sman < 0 and tman > 0: return -1

    # This reduces to direct integer comparison
    if sexp == texp: return cmp(sman, tman)

    # Check position of the highest set bit in each number. If
    # different, there is certainly an inequality.
    a = sbc + sexp
    b = tbc + texp
    if sman > 0:
        if a < b: return -1
        if a > b: return 1
    else:
        if a < b: return 1
        if a > b: return -1

    # Both numbers have the same highest bit. Subtract to find
    # how the lower bits compare.
    return cmp(fsub(s, t, 5, ROUND_FLOOR)[0], 0)


#-----------------------------------------------------------------------------
# Arithmetic: +, -, *, /, and related operations
#

def fadd(s, t, prec, rounding):
    """Add two raw mpfs and round the result to the specified precision,
    in the specified direction."""

    # We will assume below that s has the higher exponent.
    if t[1] > s[1]:
        s, t = t, s
    sman, sexp, sbc = s
    tman, texp, tbc = t

    # Check if one operand is zero. Zero always has exp = 0; if the
    # other operand has a large exponent, its mantissa will unnecessarily
    # be shifted a huge number of bits if we don't check for this case.
    if not tman: return normalize(sman, sexp, prec, rounding)
    if not sman: return normalize(tman, texp, prec, rounding)

    #------------------------------------------------------------------------
    # More generally, if one number is huge and the other is small,
    # and in particular, if their mantissas don't overlap at all at
    # the current precision level, we can avoid work.
    #       precision
    #    |            |
    #     111111111
    #  +                 222222222
    #     ------------------------
    #  #  1111111110000...
    #
    delta = (sbc + sexp) - (tbc + texp)
    if delta > prec + 5:   # an arbitrary number ~> 3
        # The result may have to be rounded up or down. So we shift s
        # and add a dummy bit outside the precision range to force
        # rounding.
        offset = min(delta + 3, prec+3)
        sman <<= offset
        if tman > 0:
            sman += 1
        else:
            sman -= 1
        return normalize(sman, sexp-offset, prec, rounding)

    #------------------------------------------------------------------------
    #  General algorithm: we set min(s.exp, t.exp) = 0, perform exact integer
    #  addition, and then round the result.
    #                   exp = 0
    #                       |
    #                       v
    #          11111111100000   <-- s.man (padded with zeros from shifting)
    #      +        222222222   <-- t.man (no shifting necessary)
    #          --------------
    #      =   11111333333333
    #
    return normalize(tman+(sman<<(sexp-texp)), texp, prec, rounding)


def fsub(s, t, prec, rounding):
    """Return the difference of two raw mpfs, s-t. This function is
    simply a wrapper of fadd that changes the sign of t."""
    man, exp, bc = t
    return fadd(s, (-man, exp, bc), prec, rounding)


def fneg(s, prec, rounding):
    """Negate a raw mpf (return -s), rounding the result to the
    specified precision."""
    return normalize(-s[0], s[1], prec, rounding)


def fneg_exact(s):
    """Negate a raw mpf (return -s), without performing any rounding."""
    return (-s[0], s[1], s[2])


def fabs(s, prec, rounding):
    """Return abs(s) of the raw mpf s, rounded to the specified
    precision."""
    man, exp, bc = s
    if man < 0:
        return normalize(-man, exp, prec, rounding)
    return normalize(man, exp, prec, rounding)


def fmul(s, t, prec, rounding):
    """Return the product of two raw mpfs, s*t, rounded to the
    specified precision."""
    sman, sexp, sbc = s
    tman, texp, tbc = t
    # This is very simple. A possible optimization would be to throw
    # away some bits when prec is much smaller than sbc+tbc
    return normalize(sman*tman, sexp+texp, prec, rounding)


def fdiv(s, t, prec, rounding):
    """Floating-point division"""
    sman, sexp, sbc = s
    tman, texp, tbc = t

    # Same strategy as for addition: if there is a remainder, perturb
    # the result a few bits outside the precision range before rounding
    extra = max(prec - sbc + tbc + 5, 5)
    quot, rem = divmod(sman<<extra, tman)
    if rem:
        quot = (quot << 5) + 1
        extra += 5
    return normalize(quot, sexp-texp-extra, prec, rounding)


def fshift_exact(s, n):
    """Quickly multiply the raw mpf s by 2**n without rounding."""
    man, exp, bc = s
    if not man:
        return s
    return man, exp+n, bc


# TODO: use directed rounding all the way through (and, account for signs?)
def fpow(s, n, prec, rounding):
    """Compute s**n, where n is an integer"""
    n = int(n)
    if n == 0: return fone
    if n == 1: return normalize(s[0], s[1], prec, rounding)
    if n == 2: return fmul(s, s, prec, rounding)
    if n == -1: return fdiv(fone, s, prec, rounding)
    if n < 0:
        return fdiv(fone, fpow(s, -n, prec+3, ROUND_FLOOR), prec, rounding)
    # Now we perform binary exponentiation. Need to estimate precision
    # to avoid rounding from temporary operations. Roughly log_2(n)
    # operations are performed.
    prec2 = prec + int(4*math.log(n, 2) + 4)
    man, exp, bc = normalize(s[0], s[1], prec2, ROUND_FLOOR)
    pm, pe, pbc = fone
    while n:
        if n & 1:
            pm, pe, pbc = normalize(pm*man, pe+exp, prec2, ROUND_FLOOR)
            n -= 1
        man, exp, bc = normalize(man*man, exp+exp, prec2, ROUND_FLOOR)
        n = n // 2
    return normalize(pm, pe, prec, rounding)


def frand(prec):
    """Return a raw mpf chosen randomly from [0, 1), with prec bits
    in the mantissa."""
    return normalize(_random.randrange(0, 1<<prec), -prec, prec, ROUND_FLOOR)
