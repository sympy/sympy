"""
Low-level functions for arbitrary-precision floating-point arithmetic.
"""

__docformat__ = 'plaintext'

import math

from bisect import bisect
from random import getrandbits

from settings import (\
    MP_BASE, MP_ZERO, MP_ONE, MP_TWO, MP_FIVE, MODE, STRICT, gmpy,
    round_floor, round_ceiling, round_down, round_up,
    round_nearest, round_fast,
    MP_BASE_TYPE,
)

# We don't pickle tuples directly for the following reasons:
#   1: pickle uses str() for ints, which is inefficient when they are large
#   2: pickle doesn't work for gmpy mpzs
# Both problems are solved by using hex()

def to_pickable(x):
    sign, man, exp, bc = x
    return sign, hex(man)[2:], exp, bc

def from_pickable(x):
    sign, man, exp, bc = x
    return (sign, MP_BASE(man, 16), exp, bc)

class ComplexResult(ValueError):
    pass

#----------------------------------------------------------------------------#
#                    Some commonly needed float values                       #
#----------------------------------------------------------------------------#

# Regular number format:
# (-1)**sign * mantissa * 2**exponent, plus bitcount of mantissa
fzero = (0, MP_ZERO, 0, 0)
fnzero = (1, MP_ZERO, 0, 0)
fone = (0, MP_ONE, 0, 1)
fnone = (1, MP_ONE, 0, 1)
ftwo = (0, MP_ONE, 1, 1)
ften = (0, MP_FIVE, 1, 3)
fhalf = (0, MP_ONE, -1, 1)

# Arbitrary encoding for special numbers: zero mantissa, nonzero exponent
fnan = (0, MP_ZERO, -123, -1)
finf = (0, MP_ZERO, -456, -2)
fninf = (1, MP_ZERO, -789, -3)

# Was 1e1000; this is broken in Python 2.4
math_float_inf = 1e300 * 1e300

#----------------------------------------------------------------------------#
#           Various utilities related to precision and bit-fiddling          #
#----------------------------------------------------------------------------#

def giant_steps(start, target):
    """Return a list of integers ~= [start, 2*start, ..., target/2,
    target] describing suitable precision steps for Newton's method."""
    L = [target]
    while L[-1] > start*2:
        L = L + [L[-1]//2 + 1]
    return L[::-1]

def giant_steps2(start, target):
    """Return a list of integers ~= [start, 3*start, ..., target/3,
    target] describing suitable precision steps for Halley's method."""
    L = [target]
    while L[-1] > start*3:
        L = L + [L[-1]//3 + 1]
    return L[::-1]

def rshift(x, n):
    """For an integer x, calculate x >> n with the fastest (floor)
    rounding. Unlike the plain Python expression (x >> n), n is
    allowed to be negative, in which case a left shift is performed."""
    if n >= 0: return x >> n
    else:      return x << (-n)

def lshift(x, n):
    """For an integer x, calculate x << n. Unlike the plain Python
    expression (x << n), n is allowed to be negative, in which case a
    right shift with default (floor) rounding is performed."""
    if n >= 0: return x << n
    else:      return x >> (-n)

def python_trailing(n):
    """Count the number of trailing zero bits in abs(n)."""
    if not n:
        return 0
    t = 0
    while not n & 1:
        n >>= 1
        t += 1
    return t

def gmpy_trailing(n):
    """Count the number of trailing zero bits in abs(n) using gmpy."""
    if n: return MP_BASE(n).scan1()
    else: return 0

# Small powers of 2
powers = [1<<_ for _ in range(300)]

def python_bitcount(n):
    """Calculate bit size of the nonnegative integer n."""
    bc = bisect(powers, n)
    if bc != 300:
        return bc
    bc = int(math.log(n, 2)) - 4
    return bc + bctable[n>>bc]

def gmpy_bitcount(n):
    """Calculate bit size of the nonnegative integer n."""
    if n: return MP_BASE(n).numdigits(2)
    else: return 0

if MODE == 'gmpy':
    bitcount = gmpy_bitcount
    trailing = gmpy_trailing
else:
    bitcount = python_bitcount
    trailing = python_trailing

# Used to avoid slow function calls as far as possible
trailtable = map(trailing, range(256))
bctable = map(bitcount, range(1024))

#----------------------------------------------------------------------------#
#                                  Rounding                                  #
#----------------------------------------------------------------------------#


# This function can be used to round a mantissa generally. However,
# we will try to do most rounding inline for efficiency.
def round_int(x, n, rnd):
    if rnd is round_nearest:
        if x >= 0:
            t = x >> (n-1)
            if t & 1 and ((t & 2) or (x & h_mask[n<300][n])):
                return (t>>1)+1
            else:
                return t>>1
        else:
            return -round_int(-x, n, rnd)
    if rnd is round_floor:
        return x >> n
    if rnd is round_ceiling:
        return -((-x) >> n)
    if rnd is round_down:
        if x >= 0:
            return x >> n
        return -((-x) >> n)
    if rnd is round_up:
        if x >= 0:
            return -((-x) >> n)
        return x >> n

# These masks are used to pick out segments of numbers to determine
# which direction to round when rounding to nearest.
class h_mask_big:
    def __getitem__(self, n):
        return (MP_ONE<<(n-1))-1

h_mask_small = [0]+[((MP_ONE<<(_-1))-1) for _ in range(1, 300)]
h_mask = [h_mask_big(), h_mask_small]

# The >> operator rounds to floor. shifts_down[rnd][sign]
# tells whether this is the right direction to use, or if the
# number should be negated before shifting
shifts_down = {round_floor:(1,0), round_ceiling:(0,1),
    round_down:(1,1), round_up:(0,0)}


#----------------------------------------------------------------------------#
#                          Normalization of raw mpfs                         #
#----------------------------------------------------------------------------#

# This function is called almost every time an mpf is created.
# It has been optimized accordingly.

def _normalize(sign, man, exp, bc, prec, rnd):
    """
    Create a raw mpf tuple with value (-1)**sign * man * 2**exp and
    normalized mantissa. The mantissa is rounded in the specified
    direction if its size exceeds the precision. Trailing zero bits
    are also stripped from the mantissa to ensure that the
    representation is canonical.

    Conditions on the input:
    * The input must represent a regular (finite) number
    * The sign bit must be 0 or 1
    * The mantissa must be positive
    * The exponent must be an integer
    * The bitcount must be exact

    If these conditions are not met, use from_man_exp, mpf_pos, or any
    of the conversion functions to create normalized raw mpf tuples.
    """
    if not man:
        return fzero
    # Cut mantissa down to size if larger than target precision
    n = bc - prec
    if n > 0:
        if rnd is round_nearest:
            t = man >> (n-1)
            if t & 1 and ((t & 2) or (man & h_mask[n<300][n])):
                man = (t>>1)+1
            else:
                man = t>>1
        elif shifts_down[rnd][sign]:
            man >>= n
        else:
            man = -((-man)>>n)
        exp += n
        bc = prec
    # Strip trailing bits
    if not man & 1:
        t = trailtable[int(man & 255)]
        if not t:
            while not man & 255:
                man >>= 8
                exp += 8
                bc -= 8
            t = trailtable[int(man & 255)]
        man >>= t
        exp += t
        bc -= t
    # Bit count can be wrong if the input mantissa was 1 less than
    # a power of 2 and got rounded up, thereby adding an extra bit.
    # With trailing bits removed, all powers of two have mantissa 1,
    # so this is easy to check for.
    if man == 1:
        bc = 1
    return sign, man, exp, bc

def _normalize1(sign, man, exp, bc, prec, rnd):
    """same as normalize, but with the added condition that
       man is odd or zero
    """
    if not man:
        return fzero
    if bc <= prec:
        return sign, man, exp, bc
    n = bc - prec
    if rnd is round_nearest:
        t = man >> (n-1)
        if t & 1 and ((t & 2) or (man & h_mask[n<300][n])):
            man = (t>>1)+1
        else:
            man = t>>1
    elif shifts_down[rnd][sign]:
        man >>= n
    else:
        man = -((-man)>>n)
    exp += n
    bc = prec
    # Strip trailing bits
    if not man & 1:
        t = trailtable[int(man & 255)]
        if not t:
            while not man & 255:
                man >>= 8
                exp += 8
                bc -= 8
            t = trailtable[int(man & 255)]
        man >>= t
        exp += t
        bc -= t
    # Bit count can be wrong if the input mantissa was 1 less than
    # a power of 2 and got rounded up, thereby adding an extra bit.
    # With trailing bits removed, all powers of two have mantissa 1,
    # so this is easy to check for.
    if man == 1:
        bc = 1
    return sign, man, exp, bc

def strict_normalize(sign, man, exp, bc, prec, rnd):
    """Additional checks on the components of an mpf. Enable tests by setting
       the environment variable MPMATH_STRICT to Y."""
    assert type(man) == MP_BASE_TYPE
    assert type(bc) in (int, long)
    assert type(exp) in (int, long)
    assert bc == bitcount(man)
    return _normalize(sign, man, exp, bc, prec, rnd)

def strict_normalize1(sign, man, exp, bc, prec, rnd):
    """Additional checks on the components of an mpf. Enable tests by setting
       the environment variable MPMATH_STRICT to Y."""
    assert type(man) == MP_BASE_TYPE
    assert type(bc) in (int, long)
    assert type(exp) in (int, long)
    assert bc == bitcount(man)
    assert (not man) or (man & 1)
    return _normalize1(sign, man, exp, bc, prec, rnd)

if STRICT:
    normalize = strict_normalize
    normalize1 = strict_normalize1
else:
    normalize = _normalize
    normalize1 = _normalize1

#----------------------------------------------------------------------------#
#                            Conversion functions                            #
#----------------------------------------------------------------------------#

def from_man_exp(man, exp, prec=None, rnd=round_fast):
    """Create raw mpf from (man, exp) pair. The mantissa may be signed.
    If no precision is specified, the mantissa is stored exactly."""
    man = MP_BASE(man)
    sign = 0
    if man < 0:
        sign = 1
        man = -man
    if man < 1024:
        bc = bctable[int(man)]
    else:
        bc = bitcount(man)
    if not prec:
        if not man:
            return fzero
        while not man & 1:
            man >>= 1
            exp += 1
            bc -= 1
        return (sign, man, exp, bc)
    return normalize(sign, man, exp, bc, prec, rnd)

int_cache = dict((n, from_man_exp(n, 0)) for n in range(-10, 257))

def from_int(n, prec=None, rnd=round_fast):
    """Create a raw mpf from an integer. If no precision is specified,
    the mantissa is stored exactly."""
    if not prec:
        if n in int_cache:
            return int_cache[n]
    return from_man_exp(n, 0, prec, rnd)

def to_man_exp(s):
    """Return (man, exp) of a raw mpf. Raise an error if inf/nan."""
    sign, man, exp, bc = s
    if (not man) and exp:
        raise ValueError("mantissa and exponent are undefined for %s" % man)
    return man, exp

def to_int(s, rnd=None):
    """Convert a raw mpf to the nearest int. Rounding is done down by
    default (same as int(float) in Python), but can be changed. If the
    input is inf/nan, an exception is raised."""
    sign, man, exp, bc = s
    if (not man) and exp:
        raise ValueError("cannot convert %s to int" % man)
    if exp >= 0:
        if sign:
            return (-man) << exp
        return man << exp
    # Make default rounding fast
    if not rnd:
        if sign:
            return -(man >> (-exp))
        else:
            return man >> (-exp)
    if sign:
        return round_int(-man, -exp, rnd)
    else:
        return round_int(man, -exp, rnd)

def mpf_ceil(s, prec, rnd=round_fast):
    """Calculate ceil of a raw mpf, and round the result in the given
    direction (not necessarily ceiling). Note: returns a raw mpf
    representing an integer, not a Python int."""
    sign, man, exp, bc = s
    if (not man) and exp:
        return s
    if exp > 0:
        return mpf_pos(s, prec, rnd)
    return from_int(to_int(s, round_ceiling), prec, rnd)

def mpf_floor(s, prec, rnd=round_fast):
    """Calculate floor of a raw mpf, and round the result in the given
    direction (not necessarily floor). Note: returns a raw mpf
    representing an integer, not a Python int."""
    sign, man, exp, bc = s
    if (not man) and exp:
        return s
    if exp > 0:
        return mpf_pos(s, prec, rnd)
    return from_int(to_int(s, round_floor), prec, rnd)

def from_float(x, prec=53, rnd=round_fast):
    """Create a raw mpf from a Python float, rounding if necessary.
    If prec >= 53, the result is guaranteed to represent exactly the
    same number as the input. If prec is not specified, use prec=53."""
    # frexp only raises an exception for nan on some platforms
    if x != x:
        return fnan
    # in Python2.5 math.frexp gives an exception for float infinity
    # in Python2.6 it returns (float infinity, 0)
    try:
        m, e = math.frexp(x)
    except:
        if x == math_float_inf: return finf
        if x == -math_float_inf: return fninf
        return fnan
    if x == math_float_inf: return finf
    if x == -math_float_inf: return fninf
    return from_man_exp(int(m*(1<<53)), e-53, prec, rnd)

def to_float(s):
    """Convert a raw mpf to a Python float. The result is exact if the
    bitcount of s is <= 53 and no underflow/overflow occurs. An
    OverflowError is raised if the number is too large to be
    represented as a regular float."""
    sign, man, exp, bc = s
    if not man:
        if s == fzero: return 0.0
        if s == finf: return 1e1000
        if s == fninf: return -1e1000
        return 1e1000/1e1000
    if sign:
        man = -man
    if bc < 100:
        return math.ldexp(man, exp)
    # Try resizing the mantissa. Overflow may still happen here.
    n = bc - 53
    m = man >> n
    return math.ldexp(m, exp + n)

def from_rational(p, q, prec, rnd=round_fast):
    """Create a raw mpf from a rational number p/q, rnd if
    necessary."""
    return mpf_div(from_int(p), from_int(q), prec, rnd)

def to_rational(s):
    """Convert a raw mpf to a rational number. Return integers (p, q)
    such that s = p/q exactly."""
    sign, man, exp, bc = s
    if sign:
        man = -man
    if bc == -1:
        raise ValueError("cannot convert %s to a rational number" % man)
    if exp >= 0:
        return man * (1<<exp), 1
    else:
        return man, 1<<(-exp)

def to_fixed(s, prec):
    """Convert a raw mpf to a fixed-point big integer"""
    sign, man, exp, bc = s
    offset = exp + prec
    if sign:
        if offset >= 0: return (-man) << offset
        else:           return (-man) >> (-offset)
    else:
        if offset >= 0: return man << offset
        else:           return man >> (-offset)


##############################################################################
##############################################################################

#----------------------------------------------------------------------------#
#                       Arithmetic operations, etc.                          #
#----------------------------------------------------------------------------#

def mpf_rand(prec):
    """Return a raw mpf chosen randomly from [0, 1), with prec bits
    in the mantissa."""
    return from_man_exp(getrandbits(prec), -prec, prec, round_floor)

def mpf_eq(s, t):
    """Test equality of two raw mpfs. This is simply tuple comparion
    unless either number is nan, in which case the result is False."""
    if not s[1] or not t[1]:
        if s == fnan or t == fnan:
            return False
    return s == t

def mpf_hash(s):
    try:
        # Try to be compatible with hash values for floats and ints
        return hash(to_float(s))
    except OverflowError:
        # We must unfortunately sacrifice compatibility with ints here. We
        # could do hash(man << exp) when the exponent is positive, but
        # this would cause unreasonable inefficiency for large numbers.
        return hash(s)

def mpf_cmp(s, t):
    """Compare the raw mpfs s and t. Return -1 if s < t, 0 if s == t,
    and 1 if s > t. (Same convention as Python's cmp() function.)"""

    # In principle, a comparison amounts to determining the sign of s-t.
    # A full subtraction is relatively slow, however, so we first try to
    # look at the components.
    ssign, sman, sexp, sbc = s
    tsign, tman, texp, tbc = t

    # Handle zeros and special numbers
    if not sman or not tman:
        if s == fzero: return -mpf_sign(t)
        if t == fzero: return mpf_sign(s)
        if s == t: return 0
        # Follow same convention as Python's cmp for float nan
        if t == fnan: return 1
        if s == finf: return 1
        return -1
    # Different sides of zero
    if ssign != tsign:
        if not ssign: return 1
        return -1
    # This reduces to direct integer comparison
    if sexp == texp:
        if ssign: return -cmp(sman, tman)
        else:     return cmp(sman, tman)
    # Check position of the highest set bit in each number. If
    # different, there is certainly an inequality.
    a = sbc + sexp
    b = tbc + texp
    if ssign:
        if a < b: return 1
        if a > b: return -1
    else:
        if a < b: return -1
        if a > b: return 1

    # Both numbers have the same highest bit. Subtract to find
    # how the lower bits compare.
    delta = mpf_sub(s, t, 5, round_floor)
    if delta[0]:
        return -1
    return 1

def mpf_lt(s, t):
    if s == fnan or t == fnan:
        return False
    return mpf_cmp(s, t) < 0

def mpf_le(s, t):
    if s == fnan or t == fnan:
        return False
    return mpf_cmp(s, t) <= 0

def mpf_gt(s, t):
    if s == fnan or t == fnan:
        return False
    return mpf_cmp(s, t) > 0

def mpf_ge(s, t):
    if s == fnan or t == fnan:
        return False
    return mpf_cmp(s, t) >= 0

def mpf_pos(s, prec, rnd=round_fast):
    """Calculate 0+s for a raw mpf (i.e., just round s to the specified
    precision)."""
    sign, man, exp, bc = s
    if (not man) and exp:
        return s
    return normalize1(sign, man, exp, bc, prec, rnd)

def mpf_neg(s, prec=None, rnd=round_fast):
    """Negate a raw mpf (return -s), rounding the result to the
    specified precision. The prec argument can be omitted to do the
    operation exactly."""
    sign, man, exp, bc = s
    if not man:
        if exp:
            if s == finf: return fninf
            if s == fninf: return finf
        return s
    if not prec:
        return (1-sign, man, exp, bc)
    return normalize1(1-sign, man, exp, bc, prec, rnd)

def mpf_abs(s, prec=None, rnd=round_fast):
    """Return abs(s) of the raw mpf s, rounded to the specified
    precision. The prec argument can be omitted to generate an
    exact result."""
    sign, man, exp, bc = s
    if (not man) and exp:
        if s == fninf:
            return finf
        return s
    if not prec:
        if sign:
            return (not sign, man, exp, bc)
        return s
    return normalize1(0, man, exp, bc, prec, rnd)

def mpf_sign(s):
    """Return -1, 0, or 1 (as a Python int, not a raw mpf) depending on
    whether s is negative, zero, or positive. (Nan is taken to give 0.)"""
    sign, man, exp, bc = s
    if not man:
        if s == finf: return 1
        if s == fninf: return -1
        return 0
    return (-1) ** sign

def mpf_add(s, t, prec, rnd=round_fast):
    if t[2] > s[2]:
        s, t = t, s
    ssign, sman, sexp, sbc = s
    tsign, tman, texp, tbc = t

    if not sman or not tman:
        if ((not sman) and sexp) or ((not tman) and texp):
            either = s, t
            if fnan in either: return fnan
            if finf in either and fninf in either: return fnan
            if finf in either: return finf
            return fninf
        # Check if one operand is zero. Zero always has exp = 0; if the
        # other operand has a huge exponent, its mantissa will unnecessarily
        # be shifted into something huge if we don't check for this case.
        if not tman: return normalize1(ssign, sman, sexp, sbc, prec, rnd)
        if not sman: return normalize1(tsign, tman, texp, tbc, prec, rnd)

    # More generally, if one number is huge and the other is small,
    # and in particular, if their mantissas don't overlap at all at
    # the current precision level, we can avoid work.
    #         precision
    #      |            |
    #       111111111
    #    +                    222222222
    #       ------------------------
    #       1111111110000... (222)
    offset = sexp - texp
    if offset > 100:
        delta = sbc + sexp - tbc - texp
        if delta > prec + 4:
            offset = min(delta, prec) + 4
            sman <<= offset
            if tsign: sman -= 1
            else:     sman += 1
            # TODO: use that bc ~= sbc+offset
            bc = bitcount(sman)
            return normalize1(ssign, sman, sexp-offset, bc, prec, rnd)
    if offset:
        if ssign == tsign:
            man = tman + (sman << offset)
            sbc += offset
            if tbc > sbc: bc = tbc - 4
            else:         bc = sbc - 4
            if bc < 4:    bc = bctable[int(man)]
            else:         bc += bctable[int(man>>bc)]
            return normalize1(ssign, man, texp, bc, prec, rnd)
        else:
            if ssign: man = tman - (sman << offset)
            else:     man = (sman << offset) - tman
            if man >= 0:
                ssign = 0
            else:
                man = -man
                ssign = 1
            bc = bitcount(man)
            return normalize1(ssign, man, texp, bc, prec, rnd)
    else:
        if ssign == tsign:
            man = tman + sman
            if tbc > sbc: bc = tbc - 4
            else:         bc = sbc - 4
            if bc < 4:    bc = bctable[int(man)]
            else:         bc += bctable[int(man>>bc)]
            return normalize(ssign, man, texp, bc, prec, rnd)
        else:
            if ssign: man = tman - sman
            else:     man = sman - tman
            if man >= 0:
                ssign = 0
            else:
                man = -man
                ssign = 1
            bc = bitcount(man)
            return normalize(ssign, man, texp, bc, prec, rnd)

def mpf_sub(s, t, prec, rnd=round_fast):
    """Return the difference of two raw mpfs, s-t. This function is
    simply a wrapper of mpf_add that changes the sign of t."""
    sign, man, exp, bc = t
    if (not man) and exp:
        return mpf_add(s, mpf_neg(t, prec, rnd), prec, rnd)
    return mpf_add(s, (1-sign, man, exp, bc), prec, rnd)

def mpf_mul(s, t, prec=0, rnd=round_fast):
    """Multiply two raw mpfs"""
    ssign, sman, sexp, sbc = s
    tsign, tman, texp, tbc = t
    sign = ssign ^ tsign
    man = sman*tman
    if not man:
        s_special = (not sman) and sexp
        t_special = (not tman) and texp
        if not s_special and not t_special:
            return fzero
        if fnan in (s, t): return fnan
        if (not tman) and texp: s, t = t, s
        if t == fzero: return fnan
        return {1:finf, -1:fninf}[mpf_sign(s) * mpf_sign(t)]
    bc = sbc + tbc - 4
    if bc < 4: bc = bctable[int(man)]
    else:      bc += bctable[int(man>>bc)]
    if prec:
        return normalize1(sign, man, sexp+texp, bc, prec, rnd)
    else:
        return (sign, man, sexp+texp, bc)

def mpf_mul_int(s, n, prec, rnd=round_fast):
    """Multiply by a Python integer."""
    sign, man, exp, bc = s
    if not man:
        return mpf_mul(s, from_int(n), prec, rnd)
    if not n:
        return fzero
    if n < 0:
        sign ^= 1
        n = -n
    man *= n
    # Generally n will be small
    if n < 1024:
        bc += bctable[int(n)] - 4
    else:
        bc += bitcount(n) - 4
    if bc < 4: bc = bctable[int(man)]
    else:      bc += bctable[int(man>>bc)]
    return normalize(sign, man, exp, bc, prec, rnd)

def mpf_shift(s, n):
    """Quickly multiply the raw mpf s by 2**n without rounding."""
    sign, man, exp, bc = s
    if not man:
        return s
    return sign, man, exp+n, bc

def mpf_frexp(x):
    """Convert x = y*2**n to (y, n) with abs(y) in [0.5, 1) if nonzero"""
    sign, man, exp, bc = x
    if not man:
        if x == fzero:
            return (fzero, 0)
        else:
            raise ValueError
    return mpf_shift(x, -bc-exp), bc+exp

def mpf_div(s, t, prec, rnd=round_fast):
    """Floating-point division"""
    ssign, sman, sexp, sbc = s
    tsign, tman, texp, tbc = t
    if not sman or not tman:
        if s == fzero:
            if t == fzero: raise ZeroDivisionError
            if t == fnan: return fnan
            return fzero
        if t == fzero:
            raise ZeroDivisionError
        s_special = (not sman) and sexp
        t_special = (not tman) and texp
        if s_special and t_special:
            return fnan
        if s == fnan or t == fnan:
            return fnan
        if not t_special:
            if t == fzero:
                return fnan
            return {1:finf, -1:fninf}[mpf_sign(s) * mpf_sign(t)]
        return fzero
    if tman == 1:
        return normalize1(ssign^tsign, sman, sexp-texp, sbc, prec, rnd)
    if ssign:
        sman = -sman
    if tsign:
        tman = -tman
    # Same strategy as for addition: if there is a remainder, perturb
    # the result a few bits outside the precision range before rounding
    extra = prec - sbc + tbc + 5
    if extra < 5:
        extra = 5
    quot, rem = divmod(sman<<extra, tman)
    if quot >= 0:
        sign = 0
    else:
        quot = -quot
        sign = 1
    if rem:
        quot = (quot << 5) + 1
        extra += 5
    bc = sbc+extra-tbc-4
    if bc < 4: bc = bctable[int(quot)]
    else:      bc += bctable[int(quot>>bc)]
    return normalize(sign, quot, sexp-texp-extra, bc, prec, rnd)

def mpf_rdiv_int(n, t, prec, rnd=round_fast):
    """Floating-point division with a Python integer as numerator"""
    tsign, tman, texp, tbc = t
    if not n or not tman:
        return mpf_div(from_int(n), t, prec, rnd)
    if tsign:
        tman = -tman
    extra = prec + tbc + 5
    quot, rem = divmod(n<<extra, tman)
    if quot >= 0:
        sign = 0
    else:
        quot = -quot
        sign = 1
    if rem:
        quot = (quot << 5) + 1
        extra += 5
        return normalize1(sign, quot, -texp-extra, bitcount(quot), prec, rnd)
    return normalize(sign, quot, -texp-extra, bitcount(quot), prec, rnd)

def mpf_mod(s, t, prec, rnd=round_fast):
    ssign, sman, sexp, sbc = s
    tsign, tman, texp, tbc = t
    if ((not sman) and sexp) or ((not tman) and texp):
        return fnan
    # Important special case: do nothing if t is larger
    if ssign == tsign and texp > sexp+sbc:
        return s
    # Another important special case: this allows us to do e.g. x % 1.0
    # to find the fractional part of x, and it will work when x is huge.
    if tman == 1 and sexp > texp+tbc:
        return fzero
    base = min(sexp, texp)
    sman = (-1)**ssign * sman
    tman = (-1)**tsign * tman
    man = (sman << (sexp-base)) % (tman << (texp-base))
    if man >= 0:
        sign = 0
    else:
        man = -man
        sign = 1
    return normalize(sign, man, base, bitcount(man), prec, rnd)

reciprocal_rnd = {
  round_down : round_up,
  round_up : round_down,
  round_floor : round_ceiling,
  round_ceiling : round_floor,
  round_nearest : round_nearest
}

negative_rnd = {
  round_down : round_down,
  round_up : round_up,
  round_floor : round_ceiling,
  round_ceiling : round_floor,
  round_nearest : round_nearest
}

def mpf_pow_int(s, n, prec, rnd=round_fast):
    """Compute s**n, where s is a raw mpf and n is a Python integer."""
    sign, man, exp, bc = s

    if (not man) and exp:
        if s == finf:
            if n > 0: return s
            if n == 0: return fnan
            return fzero
        if s == fninf:
            if n > 0: return [finf, fninf][n & 1]
            if n == 0: return fnan
            return fzero
        return fnan

    n = int(n)
    if n == 0: return fone
    if n == 1: return mpf_pos(s, prec, rnd)
    if n == 2:
        _, man, exp, bc = s
        if not man:
            return fzero
        man = man*man
        if man == 1:
            return (0, MP_ONE, exp+exp, 1)
        bc = bc + bc - 2
        bc += bctable[int(man>>bc)]
        return normalize1(0, man, exp+exp, bc, prec, rnd)
    if n == -1: return mpf_div(fone, s, prec, rnd)
    if n < 0:
        inverse = mpf_pow_int(s, -n, prec+5, reciprocal_rnd[rnd])
        return mpf_div(fone, inverse, prec, rnd)

    result_sign = sign & n

    # Use exact integer power when the exact mantissa is small
    if man == 1:
        return (result_sign, MP_ONE, exp*n, 1)
    if bc*n < 1000:
        man **= n
        return normalize1(result_sign, man, exp*n, bitcount(man), prec, rnd)

    # Use directed rounding all the way through to maintain rigorous
    # bounds for interval arithmetic
    rounds_down = (rnd is round_nearest) or \
        shifts_down[rnd][result_sign]

    # Now we perform binary exponentiation. Need to estimate precision
    # to avoid rounding errors from temporary operations. Roughly log_2(n)
    # operations are performed.
    workprec = prec + 4*bitcount(n) + 4
    _, pm, pe, pbc = fone
    while 1:
        if n & 1:
            pm = pm*man
            pe = pe+exp
            pbc += bc - 2
            pbc = pbc + bctable[int(pm >> pbc)]
            if pbc > workprec:
                if rounds_down:
                    pm = pm >> (pbc-workprec)
                else:
                    pm = -((-pm) >> (pbc-workprec))
                pe += pbc - workprec
                pbc = workprec
            n -= 1
            if not n:
                break
        man = man*man
        exp = exp+exp
        bc = bc + bc - 2
        bc = bc + bctable[int(man >> bc)]
        if bc > workprec:
            if rounds_down:
                man = man >> (bc-workprec)
            else:
                man = -((-man) >> (bc-workprec))
            exp += bc - workprec
            bc = workprec
        n = n // 2

    return normalize(result_sign, pm, pe, pbc, prec, rnd)


def mpf_perturb(x, eps_sign, prec, rnd):
    """
    For nonzero x, calculate x + eps with directed rounding, where
    eps < prec relatively and eps has the given sign (0 for
    positive, 1 for negative).

    With rounding to nearest, this is taken to simply normalize
    x to the given precision.
    """
    if rnd is round_nearest:
        return mpf_pos(x, prec, rnd)
    sign, man, exp, bc = x
    eps = (eps_sign, MP_ONE, exp+bc-prec-1, 1)
    if sign:
        away = (rnd in (round_down, round_ceiling)) ^ eps_sign
    else:
        away = (rnd in (round_up, round_ceiling)) ^ eps_sign
    if away:
        return mpf_add(x, eps, prec, rnd)
    else:
        return mpf_pos(x, prec, rnd)


##############################################################################
##############################################################################

#----------------------------------------------------------------------------#
#                              Radix conversion                              #
#----------------------------------------------------------------------------#

# TODO: speed up for bases 2, 4, 8, 16, ...

def bin_to_radix(x, xbits, base, bdigits):
    """Changes radix of a fixed-point number; i.e., converts
    x * 2**xbits to floor(x * 10**bdigits)."""
    return x * (MP_BASE(base)**bdigits) >> xbits

stddigits = '0123456789abcdefghijklmnopqrstuvwxyz'

def small_numeral(n, base=10, digits=stddigits):
    """Return the string numeral of a positive integer in an arbitrary
    base. Most efficient for small input."""
    if base == 10:
        return str(n)
    digs = []
    while n:
        n, digit = divmod(n, base)
        digs.append(digits[digit])
    return "".join(digs[::-1])

def numeral_python(n, base=10, size=0, digits=stddigits):
    """Represent the integer n as a string of digits in the given base.
    Recursive division is used to make this function about 3x faster
    than Python's str() for converting integers to decimal strings.

    The 'size' parameters specifies the number of digits in n; this
    number is only used to determine splitting points and need not be
    exact."""
    if n < 0:
        return "-" + numeral(-n, base, size, digits)
    # Fast enough to do directly
    if size < 250:
        return small_numeral(n, base, digits)
    # Divide in half
    half = (size // 2) + (size & 1)
    A, B = divmod(n, base**half)
    ad = numeral(A, base, half, digits)
    bd = numeral(B, base, half, digits).rjust(half, "0")
    return ad + bd

def numeral_gmpy(n, base=10, size=0, digits=stddigits):
    """Represent the integer n as a string of digits in the given base.
    Recursive division is used to make this function about 3x faster
    than Python's str() for converting integers to decimal strings.

    The 'size' parameters specifies the number of digits in n; this
    number is only used to determine splitting points and need not be
    exact."""
    if n < 0:
        return "-" + numeral(-n, base, size, digits)
    # gmpy.digits() may cause a segmentation fault when trying to convert
    # extremely large values to a string. The size limit may need to be
    # adjusted on some platforms, but 1500000 works on Windows and Linux.
    if size < 1500000:
        return gmpy.digits(n, base)
    # Divide in half
    half = (size // 2) + (size & 1)
    A, B = divmod(n, MP_BASE(base)**half)
    ad = numeral(A, base, half, digits)
    bd = numeral(B, base, half, digits).rjust(half, "0")
    return ad + bd

if MODE == "gmpy":
    numeral = numeral_gmpy
else:
    numeral = numeral_python

def to_digits_exp(s, dps):
    """Helper function for representing the floating-point number s as
    a decimal with dps digits. Returns (sign, string, exponent) where
    sign is '' or '-', string is the digit string, and exponent is
    the decimal exponent as an int.

    If inexact, the decimal representation is rounded toward zero."""

    # Extract sign first so it doesn't mess up the string digit count
    if s[0]:
        sign = '-'
        s = mpf_neg(s)
    else:
        sign = ''
    _sign, man, exp, bc = s

    if not man:
        return '', '0', 0

    bitprec = int(dps * math.log(10,2)) + 10

    # Cut down to size
    # TODO: account for precision when doing this
    exp_from_1 = exp + bc
    if abs(exp_from_1) > 3500:
        from libelefun import mpf_ln2, mpf_ln10
        # Set b = int(exp * log(2)/log(10))
        # If exp is huge, we must use high-precision arithmetic to
        # find the nearest power of ten
        expprec = bitcount(abs(exp)) + 5
        tmp = from_int(exp)
        tmp = mpf_mul(tmp, mpf_ln2(expprec))
        tmp = mpf_div(tmp, mpf_ln10(expprec), expprec)
        b = to_int(tmp)
        s = mpf_div(s, mpf_pow_int(ften, b, bitprec), bitprec)
        _sign, man, exp, bc = s
        exponent = b
    else:
        exponent = 0

    # First, calculate mantissa digits by converting to a binary
    # fixed-point number and then converting that number to
    # a decimal fixed-point number.
    fixprec = max(bitprec - exp - bc, 0)
    fixdps = int(fixprec / math.log(10,2) + 0.5)
    sf = to_fixed(s, fixprec)
    sd = bin_to_radix(sf, fixprec, 10, fixdps)
    digits = numeral(sd, base=10, size=dps)

    exponent += len(digits) - fixdps - 1
    return sign, digits, exponent

def to_str(s, dps, strip_zeros=True):
    """Convert a raw mpf to a decimal floating-point literal with at
    most `dps` decimal digits in the mantissa (not counting extra zeros
    that may be inserted for visual purposes).

    The literal is formatted so that it can be parsed back to a number
    by to_str, float() or Decimal()."""

    # Special numbers
    if not s[1]:
        if s == fzero: return '0.0'
        if s == finf: return '+inf'
        if s == fninf: return '-inf'
        if s == fnan: return 'nan'
        raise ValueError

    # to_digits_exp rounds to floor.
    # This sometimes kills some instances of "...00001"
    sign, digits, exponent = to_digits_exp(s, dps+3)

    if not dps:
        if digits[0] in '56789':
            exponent += 1
        digits = ".0"

    else:
        # Rounding up kills some instances of "...99999"
        if len(digits) > dps and digits[dps] in '56789' and \
            (dps < 500 or digits[dps-4:dps] == '9999'):
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
                digits = ("0"*int(-exponent)) + digits
                split = 1
            else:
                split = exponent + 1
            exponent = 0
        else:
            split = 1

        digits = (digits[:split] + "." + digits[split:])

        if strip_zeros:
            # Clean up trailing zeros
            digits = digits.rstrip('0')
            if digits[-1] == ".":
                digits += "0"

    if exponent == 0 and dps: return sign + digits
    if exponent >= 0: return sign + digits + "e+" + str(exponent)
    if exponent < 0: return sign + digits + "e" + str(exponent)

def str_to_man_exp(x, base=10):
    """Helper function for from_str."""
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
    x = MP_BASE(int(x, base))
    return x, exp

special_str = {'inf':finf, '+inf':finf, '-inf':fninf, 'nan':fnan}

def from_str(x, prec, rnd=round_fast):
    """Create a raw mpf from a decimal literal, rounding in the
    specified direction if the input number cannot be represented
    exactly as a binary floating-point number with the given number of
    bits. The literal syntax accepted is the same as for Python
    floats.

    TODO: the rounding does not work properly for large exponents.
    """
    if x in special_str:
        return special_str[x]

    man, exp = str_to_man_exp(x, base=10)

    # XXX: appropriate cutoffs & track direction
    # note no factors of 5
    if abs(exp) > 400:
        s = from_int(man, prec+10)
        s = mpf_mul(s, mpf_pow_int(ften, exp, prec+10), prec, rnd)
    else:
        if exp >= 0:
            s = from_int(man * 10**exp, prec, rnd)
        else:
            s = from_rational(man, 10**-exp, prec, rnd)
    return s

# Binary string conversion. These are currently mainly used for debugging
# and could use some improvement in the future

def from_bstr(x):
    man, exp = str_to_man_exp(x, base=2)
    man = MP_BASE(man)
    sign = 0
    if man < 0:
        man = -man
        sign = 1
    bc = bitcount(man)
    return normalize(sign, man, exp, bc, bc, round_floor)

def to_bstr(x):
    sign, man, exp, bc = x
    return ['','-'][sign] + numeral(man, size=bitcount(man), base=2) + ("e%i" % exp)





##############################################################################
##############################################################################

#----------------------------------------------------------------------------#
#                                Square roots                                #
#----------------------------------------------------------------------------#


def sqrt_initial(y, prec):
    """Given y = floor(x * 2**prec), compute floor(sqrt(x) * 2**50),
    i.e. calculate a 50-bit approximation of the square root. This is
    done quickly using regular floating-point arithmetic. It is
    assumed that x ~= 1."""

    # Two cases; the second avoids overflow
    if prec < 200: return MP_BASE(y**0.5 * 2.0**(50 - prec*0.5))
    else:          return MP_BASE((y >> (prec-100))**0.5)

# XXX: doesn't work
def invsqrt_initial(y, prec):
    """Like sqrt_initial, but computes 1/sqrt(y) instead of sqrt(y)."""
    if prec < 200: return MP_BASE(y**-0.5 * 2.0**(50 + prec*0.5))
    else:          return MP_BASE((y >> (prec-100)) ** -0.5)



# We start from a 50-bit estimate for r generated with ordinary
# floating-point arithmetic, and then refines the value to full
# accuracy using the iteration

#             1  /        y  \
#    r     = --- | r  +  --- |
#     n+1     2  \ n     r_n /

# which is simply Newton's method applied to the equation r**2 = y.

# Newton's method doubles the accuracy with each step. We make use of
# this fact by only using a precision corresponding to the current
# accuracy during intermediate iterations. For example, with a 50-bit
# accurate r_1, r_2 can be computed using 100-bit precision, r_3
# using 200-bit precision, and so on. (In practice, the precision
# levels must be chosen slightly more conservatively to account for
# rounding errors in the last one or two bits.)

# It is assumed that x ~= 1 (the main mpf_sqrt() function fiddles with
# the exponent of the input to reduce it to unit magnitude before
# passing it here.)

# TODO: it would be possible to specify separate precision levels
# for the input and output, which could be useful when calculating
# pure-integer square roots.

# This is the reference implementation. In sqrt_fixed below we
# do the steps inline at low precision for a ~15% speedup
def _sqrt_fixed(y, prec):
    r = sqrt_initial(y, prec)
    extra = 10
    prevp = 50
    for p in giant_steps(50, prec+extra):
        # In the first term, we shift by the appropriate number of bits to
        # convert r from the previous precision to the current one.
        # The "-1" divides by two in the same step.

        # In the second term, we do a fixed-point division using the usual
        # formula (y<<precision)//r. The precision term is a bit messy and
        # takes into account the fact that y, r_n and r_{n+1} all have
        # different precision levels. As before, the "-1" divides by two.
        r = lshift(r, p-prevp-1) + (lshift(y, p+prevp-prec-1)//r)
        prevp = p
    return r >> extra

def gmpy_sqrt_fixed(y, prec, shifted=True):
    if shifted:
        return gmpy.sqrt(y << prec)
    else:
        extra = 10
        return gmpy.sqrt(y << (prec + 2 * extra)), extra

def python_sqrt_fixed(y, prec, shifted=True):
    """
    Square root of a fixed-point number. Given the big integer
    y = floor(x * 2**prec), this function returns floor(r * 2**prec)
    where r = sqrt(x).

    It is assumed that x ~= 1.
    """
    r = sqrt_initial(y, prec)
    extra = 10
    prevp = 50
    prec2 = prec + extra
    # unwind giant_steps; for prec2 <= 100 giant_steps(50, prec2)
    # has one element, for 100 < prec2 < 200 it has at 2 elements
    if prec2 <= 100:
        r = lshift(r, prec2-51) + ((y << 59)//r)
    elif prec2 <= 199:
        p = prec2//2 + 1
        r = (r << (p-51)) + (lshift(y, p+59-prec2)//r)
        r = (r << (prec2-p-1)) + ((y << (p+9))//r)
    else:
        prevp1 = prec2//2 + 1
        for p in giant_steps(50, prevp1):
            r = (r << (p-prevp-1)) + (y >> (prec2-p-prevp-9))//r
            prevp = p
        r = (r << (prec2-prevp1-1)) + (y << (prevp1+9))//r
    if shifted: return r >> extra
    else:       return r, extra

def gmpy_sqrt_fixed2(y, prec):
    return gmpy.sqrt(y << prec)

def python_sqrt_fixed2(y, prec):
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
    single division at high precision levels.
    """

    # XXX
    r = to_float(from_man_exp(y, -prec, 64)) ** -0.5
    r = int(r * 2**50)
    # r = invsqrt_initial(y, prec)

    extra = int(10 + 2*prec**0.3)
    prevp = 50

    for p in giant_steps(50, prec+extra):

        # This is even messier than in sqrt_fixed. As many shifts as possible
        # have been combined together for optimal speed, at a slight expense
        # of legibility.

        # Compute r**2 at precision p.
        r2 = rshift(r*r, 2*prevp - p)

        # A = r, converted from precision prevp to p
        A = lshift(r, p-prevp)

        # S = y * r2, computed at precision p. We shift y by '-prec' to
        # account for its initial precision, and by 'p' for the fixed-point
        # multiplication
        S = (lshift(y, p-prec) * r2) >> p

        # B = (3-S) and finally the outer product, both done at precision p
        B = (3<<p) - S
        r = (A*B) >> (p+1)

        prevp = p

    # Finally, sqrt(y) = y * (1/sqrt(y))
    r = (r * y) >> prec

    return r >> extra

if MODE == 'gmpy':
    sqrt_fixed = gmpy_sqrt_fixed
    sqrt_fixed2 = gmpy_sqrt_fixed2
else:
    sqrt_fixed = python_sqrt_fixed
    sqrt_fixed2 = python_sqrt_fixed2

def mpf_sqrt(s, prec, rnd=round_fast):
    """Compute the square root of a raw mpf.

    Returns a tuple representing the square root of s, rounded to the
    nearest floating-point number in the specified rounding direction.
    The input must be a tuple representing a nonnegative floating-point
    number."""

    if s == fone:
        return fone
    sign, man, exp, bc = s
    if sign:
        raise ComplexResult("square root of a negative number")
    if not man:
        return s

    # Convert to a fixed-point number with prec2 bits. Adjust
    # exponents to be even so that they can be divided in half
    prec2 = prec + 12 + (prec & 1)

    if exp & 1:
        exp -= 1
        man <<= 1
        bc += 1

    # Mantissa may have more bits than we need. Trim it down.
    shift = bc - prec2
    shift -= shift & 1
    man = rshift(man, shift)

    rnd_shift = 0
    if rnd == 'd' or rnd == 'f':
        rnd_shift = 1
    if prec < 20000:
        man, extra = sqrt_fixed(man+rnd_shift, prec2, False)
    else:
        man = sqrt_fixed2(man+rnd_shift, prec2)
        extra = 0

    return from_man_exp(man, ((exp+shift-prec2)>>1) - extra, prec, rnd)

def mpf_hypot(x, y, prec, rnd=round_fast):
    """Compute the Euclidean norm sqrt(x**2 + y**2) of two raw mpfs
    x and y."""
    if y == fzero: return mpf_abs(x, prec, rnd)
    if x == fzero: return mpf_abs(y, prec, rnd)
    hypot2 = mpf_add(mpf_mul(x,x), mpf_mul(y,y), prec+4)
    return mpf_sqrt(hypot2, prec, rnd)
