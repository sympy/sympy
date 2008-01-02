"""
Functions for converting between raw mpmath floating-point numbers and
other types and number representations (int, float, str...).
"""

from util import *
from floatop import *
from constants import flog2, flog10
import math


#----------------------------------------------------------------------
# Strings
#

def to_digits_exp(s, dps):
    """
    Helper function for representing the floating-point number s as a
    decimal with dps places. Returns (sign, string, exponent)
    containing '' or '-', the decimal digits as a string, and an
    integer for the decimal exponent.

    If inexact, the decimal representation is rounded toward zero.
    """

    # Extract sign so it doesn't mess up the string digit count
    sign = ''
    if s[0] < 0:
        sign = '-'
    s = fabs(s, s[2], ROUND_FLOOR)
    man, exp, bc = s

    if not man:
        return '', '0', 0

    bitprec = int(dps * math.log(10,2)) + 10

    # Cut down to size
    # TODO: account for precision when doing this
    exp_from_1 = exp + bc
    if abs(exp) > 500:
        # Set b = int(exp * log(2)/log(10))
        # If exp is huge, we must use high-precision arithmetic to
        # find the nearest power of ten
        expprec = bitcount(exp) + 5
        RF = ROUND_FLOOR
        tmp = from_int(exp, expprec, RF)
        tmp = fmul(tmp, flog2(expprec, RF), expprec, RF)
        tmp = fdiv(tmp, flog10(expprec, RF), expprec, RF)
        b = to_int(tmp)
        s = fdiv(s, fpow(ften, b, bitprec, RF), bitprec, RF)
        man, exp, bc = s
        exponent = b
    else:
        exponent = 0

    # First, calculate mantissa digits by converting to a binary
    # fixed-point number and then converting that number to
    # a decimal fixed-point number.
    fixprec = max(bitprec - exp, 0)
    fixdps = int(fixprec / math.log(10,2) + 0.5)
    sf = make_fixed(s, fixprec)
    sd = bin_to_radix(sf, fixprec, 10, fixdps)
    digits = numeral(sd, base=10, size=dps)

    exponent += len(digits) - fixdps - 1
    return sign, digits, exponent


def to_str(s, dps):
    """Convert a raw mpf to a decimal floating-point literal with at
    most `dps` decimal digits in the mantissa (not counting extra zeros
    that may be inserted for visual purposes).

    The literal is formatted so that it can be parsed back to a number
    by to_str, float() or Decimal()."""

    # to_digits_exp rounds to floor.
    # This sometimes kills some instances of "...00001"
    sign, digits, exponent = to_digits_exp(s, dps+3)

    # Rounding up kills some instances of "...99999"
    if len(digits) > dps and digits[dps] in '56789':
        digits2 = str(int(digits[:dps]) + 1)
        if len(digits2) > dps:
            digits2 = digits2[:dps]
            exponent += 1
        digits = digits2
    else:
        digits = digits[:dps]

    # Prettify numbers close to unit magnitude
    if -(dps//3) < exponent < dps:
        if exponent < 0:
            digits = ("0"*(-exponent)) + digits
            split = 1
        else:
            split = exponent + 1
        exponent = 0
    else:
        split = 1

    digits = (digits[:split] + "." + digits[split:])

    # Clean up trailing zeros
    digits = digits.rstrip('0')
    if digits[-1] == ".":
        digits += "0"

    if exponent == 0: return sign + digits
    if exponent > 0: return sign + digits + "e+" + str(exponent)
    if exponent < 0: return sign + digits + "e" + str(exponent)


def str_to_man_exp(x, base=10):
    # Verify that the input is a valid float literal
    float(x)
    # Split into mantissa, exponent
    x = x.lower()
    parts = x.split('e')
    if len(parts) == 1:
        exp = 0
    else: # == 2
        x = parts[0]
        exp = int(parts[1])
    # Look for radix point in mantissa
    parts = x.split('.')
    if len(parts) == 2:
        a, b = parts[0], parts[1].rstrip('0')
        exp -= len(b)
        x = a + b
    x = int(x, base)
    return x, exp

def from_str(x, prec, rounding):
    """Create a raw mpf from a decimal literal, rounding in the
    specified direction if the input number cannot be represented
    exactly as a binary floating-point number with the given number of
    bits. The literal syntax accepted is the same as for Python
    floats.

    TODO: the rounding does not work properly for large exponents.
    """
    man, exp = str_to_man_exp(x, base=10)

    # XXX: appropriate cutoffs & track direction
    # note no factors of 5
    if abs(exp) > 400:
        s = from_int(man, prec+10, ROUND_FLOOR)
        s = fmul(s, fpow(ften, exp, prec+10, ROUND_FLOOR), prec, rounding)
    else:
        if exp >= 0:
            s = from_int(man * 10**exp, prec, rounding)
        else:
            s = from_rational(man, 10**-exp, prec, rounding)
    return s

# Binary string conversion. These are currently mainly used for debugging
# and could use some improvement in the future

def from_bstr(x):
    man, exp = str_to_man_exp(x, base=2)
    return normalize(man, exp, bitcount(man), ROUND_FLOOR)

def to_bstr(x):
    man, exp, bc = x
    return numeral(man, size=bitcount(man), base=2) + ("e%i" % exp)


#----------------------------------------------------------------------
# Integers
#

def from_int(n, prec, rounding):
    """Create a raw mpf from an integer, rounding if necessary."""
    return normalize(n, 0, prec, rounding)

def from_int_exact(n):
    """Create a raw mpf from an integer, automatically choosing the
    precision high enough to represent the integer exactly."""
    t = trailing_zeros(n)
    n >>= t
    return (n, t, bitcount(n))

def to_int(s):
    """Convert a raw mpf to the nearest int, rounding down."""
    man, exp, bc = s
    return rshift(man, -exp, ROUND_DOWN)


#----------------------------------------------------------------------
# Regular python floats
#

def from_float(x, prec, rounding):
    """Create a raw mpf from a Python float, rounding if necessary.
    If prec >= 53, the result is guaranteed to represent exactly the
    same number as the input."""
    m, e = math.frexp(x)
    return normalize(int(m*(1<<53)), e-53, prec, rounding)

def to_float(s):
    """Convert a raw mpf to a Python float. The result is exact if the
    bitcount of s is <= 53 and no underflow/overflow occurs. An
    OverflowError is raised if the number is too large to be
    represented as a regular float."""
    man, exp, bc = s
    try:
        return math.ldexp(man, exp)
    except OverflowError:
        # Try resizing the mantissa. Overflow may still happen here.
        n = bc - 53
        m = man >> n
        return math.ldexp(m, exp + n)


#----------------------------------------------------------------------
# Rational numbers (for use by other libraries)
#

def from_rational(p, q, prec, rounding):
    """Create a raw mpf from a rational number p/q, rounding if
    necessary."""
    return fdiv(from_int_exact(p), from_int_exact(q), prec, rounding)

def to_rational(s):
    """Convert a raw mpf to a rational number. Return integers (p, q)
    such that s = p/q exactly. p and q are not reduced to lowest terms."""
    man, exp, bc = s
    if exp > 0:
        return man * 2**exp, 1
    else:
        return man, 2**-exp
