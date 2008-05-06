"""
Low-level functions for arbitrary-precision floating-point arithmetic.
"""

__docformat__ = 'plaintext'

import math
from bisect import bisect
from random import randrange

class ComplexResult(ValueError):
    pass

#----------------------------------------------------------------------------#
#                    Some commonly needed float values                       #
#----------------------------------------------------------------------------#

# Regular number format:
# (-1)**sign * mantissa * 2**exponent, plus bitcount of mantissa
fzero = (0, 0, 0, 0)
fnzero = (1, 0, 0, 0)
fone = (0, 1, 0, 1)
fnone = (1, 1, 0, 1)
ftwo = (0, 1, 1, 1)
ften = (0, 5, 1, 3)
fhalf = (0, 1, -1, 1)

# Arbitrary encoding for special numbers: zero mantissa, nonzero exponent
fnan = (0, 0, -123, -1)
finf = (0, 0, -456, -2)
fninf = (1, 0, -789, -3)


#----------------------------------------------------------------------------#
#           Various utilities related to precision and bit-fiddling          #
#----------------------------------------------------------------------------#

def prec_to_dps(n):
    """Return number of accurate decimals that can be represented
    with a precision of n bits."""
    return max(1, int(round(int(n)/3.3219280948873626)-1))

def dps_to_prec(n):
    """Return the number of bits required to represent n decimals
    accurately."""
    return max(1, int(round((int(n)+1)*3.3219280948873626)))

def repr_dps(n):
    """Return the number of decimal digits required to represent
    a number with n-bit precision so that it can be uniquely
    reconstructed from the representation."""
    dps = prec_to_dps(n)
    if dps == 15:
        return 17
    return dps + 3

def giant_steps(start, target):
    """Return a list of integers ~= [start, 2*start, ..., target/2,
    target] describing suitable precision steps for Newton's method."""
    L = [target]
    while L[-1] > start*2:
        L = L + [L[-1]//2 + 1]
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

def trailing(n):
    """Count the number of trailing zero bits in abs(n)."""
    if not n:
        return 0
    t = 0
    while not n & 1:
        n >>= 1
        t += 1
    return t

# Small powers of 2
powers = [1<<_ for _ in range(300)]

def bitcount(n):
    """Calculate bit size of the nonnegative integer n."""
    bc = bisect(powers, n)
    if bc != 300:
        return bc
    bc = int(math.log(n, 2)) - 4
    return bc + bctable[n>>bc]

# Used to avoid slow function calls as far as possible
trailtable = map(trailing, range(256))
bctable = map(bitcount, range(1024))

#----------------------------------------------------------------------------#
#                                  Rounding                                  #
#----------------------------------------------------------------------------#

# All supported rounding modes
round_nearest = intern('n')
round_floor = intern('f')
round_ceiling = intern('c')
round_up = intern('u')
round_down = intern('d')

round_fast = round_down

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
        return (1<<(n-1))-1

h_mask_small = [0]+[((1<<(_-1))-1) for _ in range(1, 300)]
h_mask = [h_mask_big(), h_mask_small]

# The >> operator rounds to floor. shifts_down[rnd][sign]
# tells whether this is the right direction to use, or if the
# number should be negated before shifting
shifts_down = {'f':(1,0), 'c':(0,1), 'd':(1,1), 'u':(0,0)}


#----------------------------------------------------------------------------#
#                          Normalization of raw mpfs                         #
#----------------------------------------------------------------------------#

# This function is called almost every time an mpf is created.
# It has been optimized accordingly.

def normalize(sign, man, exp, bc, prec, rnd):
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
    * The mxponent must be an integer
    * The bitcount must be exact

    If these conditions are not met, use from_man_exp, fpos, or any
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
        t = trailtable[man & 255]
        if not t:
            while not man & 255:
                man >>= 8
                exp += 8
                bc -= 8
            t = trailtable[man & 255]
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


#----------------------------------------------------------------------------#
#                            Conversion functions                            #
#----------------------------------------------------------------------------#

def from_man_exp(man, exp, prec=None, rnd=round_fast):
    """Create raw mpf from (man, exp) pair. The mantissa may be signed.
    If no precision is specified, the mantissa is stored exactly."""
    sign = 0
    if man < 0:
        sign = 1
        man = -man
    if man < 1024:
        bc = bctable[man]
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

def from_int(n, prec=None, rnd=round_fast):
    """Create a raw mpf from an integer. If no precision is specified,
    the mantissa is stored exactly."""
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

def fceil(s, prec, rnd=round_fast):
    """Calculate ceil of a raw mpf, and round the result in the given
    direction (not necessarily ceiling). Note: returns a raw mpf
    representing an integer, not a Python int."""
    sign, man, exp, bc = s
    if (not man) and exp:
        return s
    if exp > 0:
        return fpos(s, prec, rnd)
    return from_int(to_int(s, round_ceiling), prec, rnd)

def ffloor(s, prec, rnd=round_fast):
    """Calculate floor of a raw mpf, and round the result in the given
    direction (not necessarily floor). Note: returns a raw mpf
    representing an integer, not a Python int."""
    sign, man, exp, bc = s
    if (not man) and exp:
        return s
    if exp > 0:
        return fpos(s, prec, rnd)
    return from_int(to_int(s, round_floor), prec, rnd)

def from_float(x, prec=53, rnd=round_fast):
    """Create a raw mpf from a Python float, rounding if necessary.
    If prec >= 53, the result is guaranteed to represent exactly the
    same number as the input. If prec is not specified, use prec=53."""
    # frexp only raises an exception for nan on some platforms
    if x != x:
        return fnan
    try:
        m, e = math.frexp(x)
    except:
        if x == 1e1000: return finf
        if x == -1e1000: return fninf
        return fnan
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
    return fdiv(from_int(p), from_int(q), prec, rnd)

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

def frand(prec):
    """Return a raw mpf chosen randomly from [0, 1), with prec bits
    in the mantissa."""
    return from_man_exp(randrange(0, 1<<prec), -prec, prec, round_floor)

def feq(s, t):
    """Test equality of two raw mpfs. This is simply tuple comparion
    unless either number is nan, in which case the result is False."""
    if not s[1] or not t[1]:
        if s == fnan or t == fnan:
            return False
    return s == t

def fhash(s):
    try:
        # Try to be compatible with hash values for floats and ints
        return hash(to_float(s))
    except OverflowError:
        # We must unfortunately sacrifice compatibility with ints here. We
        # could do hash(man << exp) when the exponent is positive, but
        # this would cause unreasonable inefficiency for large numbers.
        return hash(s)

def fcmp(s, t):
    """Compare the raw mpfs s and t. Return -1 if s < t, 0 if s == t,
    and 1 if s > t. (Same convention as Python's cmp() function.)"""

    # In principle, a comparison amounts to determining the sign of s-t.
    # A full subtraction is relatively slow, however, so we first try to
    # look at the components.
    ssign, sman, sexp, sbc = s
    tsign, tman, texp, tbc = t

    # Handle zeros and special numbers
    if not sman or not tman:
        if s == fzero: return -fsign(t)
        if t == fzero: return fsign(s)
        if s is t: return 0
        # Follow same convention as Python's cmp for float nan
        if t is fnan: return 1
        if s is finf: return 1
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
    delta = fsub(s, t, 5, round_floor)
    if delta[0]:
        return -1
    return 1

def flt(s, t):
    if s == fnan or t == fnan:
        return False
    return fcmp(s, t) < 0

def fle(s, t):
    if s == fnan or t == fnan:
        return False
    return fcmp(s, t) <= 0

def fgt(s, t):
    if s == fnan or t == fnan:
        return False
    return fcmp(s, t) > 0

def fge(s, t):
    if s == fnan or t == fnan:
        return False
    return fcmp(s, t) >= 0

def fpos(s, prec, rnd=round_fast):
    """Calculate 0+s for a raw mpf (i.e., just round s to the specified
    precision)."""
    sign, man, exp, bc = s
    if (not man) and exp:
        return s
    return normalize(sign, man, exp, bc, prec, rnd)

def fneg(s, prec=None, rnd=round_fast):
    """Negate a raw mpf (return -s), rounding the result to the
    specified precision. The prec argument can be omitted to do the
    operation exactly."""
    sign, man, exp, bc = s
    if not man:
        if exp:
            if s is finf: return fninf
            if s is fninf: return finf
        return s
    if not prec:
        return (1-sign, man, exp, bc)
    return normalize(1-sign, man, exp, bc, prec, rnd)

def fabs(s, prec=None, rnd=round_fast):
    """Return abs(s) of the raw mpf s, rounded to the specified
    precision. The prec argument can be omitted to generate an
    exact result."""
    sign, man, exp, bc = s
    if (not man) and exp:
        if s is fninf:
            return finf
        return s
    if not prec:
        if sign:
            return (not sign, man, exp, bc)
        return s
    return normalize(0, man, exp, bc, prec, rnd)

def fsign(s):
    """Return -1, 0, or 1 (as a Python int, not a raw mpf) depending on
    whether s is negative, zero, or positive. (Nan is taken to give 0.)"""
    sign, man, exp, bc = s
    if not man:
        if s is finf: return 1
        if s is fninf: return -1
        return 0
    return (-1) ** sign

def fadd(s, t, prec, rnd=round_fast):
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
        if not tman: return normalize(ssign, sman, sexp, sbc, prec, rnd)
        if not sman: return normalize(tsign, tman, texp, tbc, prec, rnd)

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
            return normalize(ssign, sman, sexp-offset, bc, prec, rnd)
    if ssign == tsign:
        man = tman + (sman << offset)
        sbc += offset
        if tbc > sbc: bc = tbc - 4
        else:         bc = sbc - 4
        if bc < 4:    bc = bctable[man]
        else:         bc += bctable[man>>bc]
        return normalize(ssign, man, texp, bc, prec, rnd)
    else:
        if ssign: man = tman - (sman << offset)
        else:     man = (sman << offset) - tman
        if man >= 0:
            sign = 0
        else:
            man = -man
            sign = 1
        bc = bitcount(man)
        return normalize(sign, man, texp, bc, prec, rnd)

def fsub(s, t, prec, rnd=round_fast):
    """Return the difference of two raw mpfs, s-t. This function is
    simply a wrapper of fadd that changes the sign of t."""
    sign, man, exp, bc = t
    if (not man) and exp:
        return fadd(s, fneg(t, prec, rnd), prec, rnd)
    return fadd(s, (1-sign, man, exp, bc), prec, rnd)

def fmul(s, t, prec=0, rnd=round_fast):
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
        return {1:finf, -1:fninf}[fsign(s) * fsign(t)]
    bc = sbc + tbc - 4
    if bc < 4: bc = bctable[man]
    else:      bc += bctable[man>>bc]
    if prec:
        return normalize(sign, man, sexp+texp, bc, prec, rnd)
    else:
        return (sign, man, sexp+texp, bc)

def fmuli(s, n, prec, rnd=round_fast):
    """Multiply by a Python integer."""
    sign, man, exp, bc = s
    if not man:
        return fmul(s, from_int(n), prec, rnd)
    if not n:
        return fzero
    if n < 0:
        sign ^= 1
        n = -n
    man *= n
    # Generally n will be small
    if n < 1024:
        bc += bctable[n] - 4
    else:
        bc += bitcount(n) - 4
    if bc < 4: bc = bctable[man]
    else:      bc += bctable[man>>bc]
    return normalize(sign, man, exp, bc, prec, rnd)

def fshift(s, n):
    """Quickly multiply the raw mpf s by 2**n without rnd."""
    sign, man, exp, bc = s
    if not man:
        return s
    return sign, man, exp+n, bc

def fdiv(s, t, prec, rnd=round_fast):
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
            return {1:finf, -1:fninf}[fsign(s) * fsign(t)]
        return fzero
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
    if bc < 4: bc = bctable[quot]
    else:      bc += bctable[quot>>bc]
    return normalize(sign, quot, sexp-texp-extra, bc, prec, rnd)

def fdivi(n, t, prec, rnd=round_fast):
    """Floating-point division with a Python integer as numerator"""
    tsign, tman, texp, tbc = t
    if not n or not tman:
        return fdiv(from_int(n), t, prec, rnd)
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
    return normalize(sign, quot, -texp-extra, bitcount(quot), prec, rnd)

def fmod(s, t, prec, rnd=round_fast):
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

def fpow(s, t, prec, rnd=round_fast):
    """Compute s**t. Raises ComplexResult if s is negative and t is
    fractional."""
    ssign, sman, sexp, sbc = s
    tsign, tman, texp, tbc = t
    if ssign and texp < 0:
        raise ComplexResult("negative number raised to a fractional power")
    if texp >= 0:
        return fpowi(s, (-1)**tsign * (tman<<texp), prec, rnd)
    # s**(n/2) = sqrt(s)**n
    if texp == -1:
        if tman == 1:
            if tsign:
                return fdiv(fone, fsqrt(s, prec+10, reciprocal_rnd[rnd]), prec, rnd)
            return fsqrt(s, prec, rnd)
        else:
            if tsign:
                return fpowi(fsqrt(s, prec+10, reciprocal_rnd[rnd]), -tman, prec, rnd)
            return fpowi(fsqrt(s, prec+10, rnd), tman, prec, rnd)
    # General formula: s**t = exp(t*log(s))
    # TODO: handle rnd direction of the logarithm carefully
    c = flog(s, prec+10, rnd)
    return fexp(fmul(t, c), prec, rnd)

def fpowi(s, n, prec, rnd=round_fast):
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
    if n == 1: return fpos(s, prec, rnd)
    if n == 2:
        _, man, exp, bc = s
        if not man:
            return fzero
        man = man*man
        if man == 1:
            return (0, 1, exp+exp, 1)
        bc = bc + bc - 2
        bc += bctable[man>>bc]
        return normalize(0, man, exp+exp, bc, prec, rnd)
    if n == -1: return fdiv(fone, s, prec, rnd)
    if n < 0:
        inverse = fpowi(s, -n, prec+5, reciprocal_rnd[rnd])
        return fdiv(fone, inverse, prec, rnd)

    result_sign = sign & n

    # Use exact integer power when the exact mantissa is small
    if man == 1:
        return (result_sign, 1, exp*n, 1)
    if bc*n < 1000:
        man **= n
        return normalize(result_sign, man, exp*n, bitcount(man), prec, rnd)

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
            pbc = pbc + bctable[pm >> pbc]
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
        bc = bc + bctable[man >> bc]
        if bc > workprec:
            if rounds_down:
                man = man >> (bc-workprec)
            else:
                man = -((-man) >> (bc-workprec))
            exp += bc - workprec
            bc = workprec
        n = n // 2

    return normalize(result_sign, pm, pe, pbc, prec, rnd)


##############################################################################
##############################################################################

#----------------------------------------------------------------------------#
#                              Radix conversion                              #
#----------------------------------------------------------------------------#

# TODO: speed up for bases 2, 4, 8, 16, ...

def bin_to_radix(x, xbits, base, bdigits):
    """Changes radix of a fixed-point number; i.e., converts
    x * 2**xbits to floor(x * 10**bdigits)."""
    return x * (base**bdigits) >> xbits

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

def numeral(n, base=10, size=0, digits=stddigits):
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

def to_digits_exp(s, dps):
    """Helper function for representing the floating-point number s as
    a decimal with dps digits. Returns (sign, string, exponent) where
    sign is '' or '-', string is the digit string, and exponent is
    the decimal exponent as an int.

    If inexact, the decimal representation is rounded toward zero."""

    # Extract sign first so it doesn't mess up the string digit count
    if s[0]:
        sign = '-'
        s = fneg(s)
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
        # Set b = int(exp * log(2)/log(10))
        # If exp is huge, we must use high-precision arithmetic to
        # find the nearest power of ten
        expprec = bitcount(abs(exp)) + 5
        tmp = from_int(exp)
        tmp = fmul(tmp, flog2(expprec))
        tmp = fdiv(tmp, flog10(expprec), expprec)
        b = to_int(tmp)
        s = fdiv(s, fpowi(ften, b, bitprec), bitprec)
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

def to_str(s, dps):
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
    x = int(x, base)
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
        s = fmul(s, fpowi(ften, exp, prec+10), prec, rnd)
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
    if prec < 200: return int(y**0.5 * 2.0**(50 - prec*0.5))
    else:          return int((y >> (prec-100))**0.5)

# XXX: doesn't work
def invsqrt_initial(y, prec):
    """Like sqrt_initial, but computes 1/sqrt(y) instead of sqrt(y)."""
    if prec < 200: return int(y**-0.5 * 2.0**(50 + prec*0.5))
    else:          return int((y >> (prec-100)) ** -0.5)



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

# It is assumed that x ~= 1 (the main fsqrt() function fiddles with
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

def sqrt_fixed(y, prec,shifted=True):
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

def fsqrt(s, prec, rnd=round_fast):
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
        return fzero

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

    if prec < 20000:
        man, extra = sqrt_fixed(man, prec2, False)
    else:
        man = sqrt_fixed2(man, prec2)
        extra = 0

    return from_man_exp(man, ((exp+shift-prec2)>>1) - extra, prec, rnd)

def fhypot(x, y, prec, rnd=round_fast):
    """Compute the Euclidean norm sqrt(x**2 + y**2) of two raw mpfs
    x and y."""
    if y == fzero: return fabs(x, prec, rnd)
    if x == fzero: return fabs(y, prec, rnd)
    hypot2 = fadd(fmul(x,x), fmul(y,y), prec+4)
    return fsqrt(hypot2, prec, rnd)


##############################################################################
##############################################################################

#----------------------------------------------------------------------------#
#                           Mathematical constants                           #
#----------------------------------------------------------------------------#

def constant_memo(f):
    """Cache computed values of mathematical constants"""
    f.memo_prec = -1
    f.memo_val = None
    def g(prec):
        if prec == f.memo_prec:
            return f.memo_val
        if prec < f.memo_prec:
            return f.memo_val >> (f.memo_prec-prec)
        f.memo_val = f(prec)
        f.memo_prec = prec
        return f.memo_val
    g.__name__ = f.__name__
    g.__doc__ = f.__doc__
    return g

def acot(n, prec, hyperbolic):
    """Compute acot of an integer using fixed-point arithmetic. With
    hyperbolic=True, compute acoth. The standard Taylor series
    is used."""
    s = t = (1 << prec) // n  # 1 / n
    k = 3
    while 1:
        # Repeatedly divide by k * n**2, and add
        t //= (n*n)
        term = t // k
        if not term:
            break
        # Alternate signs
        if hyperbolic or not k & 2:
            s += term
        else:
            s -= term
        k += 2
    return s

def machin(coefs, prec, hyperbolic=False):
    """Evaluate a Machin-like formula, i.e., a linear combination of
    acot(n) or acoth(n) for specific integer values of n, using fixed-
    point arithmetic. The input should be a list [(c, n), ...], giving
    c*acot[h](n) + ..."""
    extraprec = 10
    s = 0
    for a, b in coefs:
        s += a * acot(b, prec+extraprec, hyperbolic)
    return (s >> extraprec)

#----------------------------------------------------------------------------
# Pi

def agm_status(prec, step, adiff, verbose_base):
    logdiff = math.log(max(1, adiff), verbose_base)
    digits = int(prec/math.log(verbose_base,2) - logdiff)
    print "  iteration", step, ("(accuracy ~= %i base-%i digits)" % \
       (digits, verbose_base))

def pi_agm(prec, verbose=False, verbose_base=10):
    """
    Compute floor(pi * 2**prec) as a big integer using the Brent-
    Salamin algorithm based on the arithmetic-geometric mean.

    See for example Wikipedia (http://en.wikipedia.org/wiki/Brent-
    Salamin_algorithm) or "Pi and the AGM" by Jonathan and Peter
    Borwein (Wiley, 1987). The algorithm (as stated in the Wikipedia
    article) consists of setting

      a_0 = 1
      b_0 = 1/sqrt(2)
      t_0 = 1/4
      p_0 = 1

    and computing

      a_{n+1} = (a_n + b_n)/2
      b_{n+1} = sqrt(a_n * b_n)
      t_{n+1} = t_n - p_n*(a_n - a_{n+1})**2
      p_{n+1} = 2*p_n

    for n = 0, 1, 2, 3, ..., after which the approximation is given by
    pi ~= (a_n + b_n)**2 / (4*t_n). Each step roughly doubles the
    number of correct digits.
    """
    extraprec = 50
    prec += extraprec

    # Initialial values. a, b and t are fixed-point numbers
    a = 1 << prec
    b = sqrt_fixed2(a >> 1, prec)
    t = a >> 2
    p = 1

    step = 1
    while 1:
        an = (a + b) >> 1
        adiff = a - an
        if verbose:
            agm_status(prec, step, adiff, verbose_base)
        # No change in a
        if p > 16 and abs(adiff) < 1000:
            break
        prod = (a * b) >> prec
        b = sqrt_fixed2(prod, prec)
        t = t - p*((adiff**2) >> prec)
        p = 2*p
        a = an
        step += 1
    if verbose:
        print "  final division"
    pi = ((((a+b)**2) >> 2) // t)
    return pi >> extraprec

@constant_memo
def pi_fixed(prec):
    """
    Compute floor(pi * 2**prec) as a big integer.

    For low precisions, Machin's formula pi = 16*acot(5)-4*acot(239)
    is used. For high precisions, the more efficient arithmetic-
    geometric mean iteration is used.
    """
    if prec < 2000:
        return machin([(16, 5), (-4, 239)], prec)
    return pi_agm(prec)

def fpi(prec, rnd=round_fast):
    """Compute a floating-point approximation of pi"""
    return from_man_exp(pi_fixed(prec+5), -prec-5, prec, rnd)

def fdegree(prec, rnd=round_fast):
    """Compute 1 degree = pi / 180."""
    return from_man_exp(pi_fixed(prec+10)//180, -prec-10, prec, rnd)

#----------------------------------------------------------------------------
# Logarithms of integers are needed for various computations involving
# logarithms, powers, radix conversion, etc
@constant_memo
def log2_fixed(prec):
    return machin([(18, 26), (-2, 4801), (8, 8749)], prec, True)

def flog2(prec, rnd=round_fast):
    return from_man_exp(log2_fixed(prec+5), -prec-5, prec, rnd)

@constant_memo
def log10_fixed(prec):
    return machin([(46, 31), (34, 49), (20, 161)], prec, True)

def flog10(prec, rnd=round_fast):
    return from_man_exp(log10_fixed(prec+5), -prec-5, prec, rnd)

#----------------------------------------------------------------------------
# exp(1) is computed using the Taylor series for exp
@constant_memo
def e_fixed(prec):
    a = 1 << prec
    s = a << 1
    n = 2
    while a:
        a //= n
        s += a
        n += 1
    return s

def fe(prec, rnd=round_fast):
    return from_man_exp(e_fixed(prec+15), -prec-15, prec, rnd)


##############################################################################
##############################################################################

#----------------------------------------------------------------------------#
#                                                                            #
#                           Exponential function                             #
#                                                                            #
#----------------------------------------------------------------------------#

# The exponential function has a rapidly convergent Maclaurin series:
# 
#     exp(x) = 1 + x + x**2/2! + x**3/3! + x**4/4! + ...
# 
# The series can be summed very easily using fixed-point arithmetic.
# The convergence can be improved further, using a trick due to
# Richard P. Brent: instead of computing exp(x) directly, we choose a
# small integer r (say, r=10) and compute exp(x/2**r)**(2**r).

# The optimal value for r depends on the Python platform, the magnitude
# of x and the target precision, and has to be estimated from
# experimental timings. One test with x ~= 0.3 showed that
# r = 2.2*prec**0.42 gave a good fit to the optimal values for r for
# prec between 1 and 10000 bits, on one particular machine.

# This optimization makes the summation about twice as fast at
# low precision levels and much faster at high precision
# (roughly five times faster at 1000 decimal digits).

# If |x| is very large, we first rewrite it as t + n*log(2) with the
# integer n chosen such that |t| <= log(2), and then calculate
# exp(x) as exp(t)*(2**n), using the Maclaurin series for exp(t)
# (the multiplication by 2**n just amounts to shifting the exponent).

# Input: x * 2**prec
# Output: exp(x) * 2**(prec + r)
def exp_series(x, prec, r):
    x >>= r
    s = (1 << prec) + x
    a = x
    k = 2
    # Sum exp(x/2**r)
    while 1:
        a = ((a*x) >> prec) // k
        if not a:
            break
        s += a
        k += 1
    # Calculate s**(2**r) by repeated squaring
    while r:
        s = (s*s) >> prec
        r -= 1
    return s

def fexp(x, prec, rnd=round_fast):
    sign, man, exp, bc = x
    if not man:
        if not exp:
            return fone
        if x == fninf:
            return fzero
        return x
    # Fast handling e**n. TODO: the best cutoff depends on both the
    # size of n and the precision.
    if prec > 600 and exp >= 0:
        return fpowi(fe(prec+10), (-1)**sign *(man<<exp), prec, rnd)
    # extra precision needs to be similar in magnitude to log_2(|x|)
    # for the modulo reduction, plus r for the error from squaring r times
    wp = prec + max(0, bc+exp)
    r = int(2 * wp**0.4)
    wp += r + 20
    t = to_fixed(x, wp)
    # abs(x) > 1?
    if exp+bc > 1:
        lg2 = log2_fixed(wp)
        n, t = divmod(t, lg2)
    else:
        n = 0
    man = exp_series(t, wp, r)
    bc = wp - 2 + bctable[man >> (wp - 2)]
    return normalize(0, man, -wp+n, bc, prec, rnd)


#----------------------------------------------------------------------------#
#                                                                            #
#                                Logarithms                                  #
#                                                                            #
#----------------------------------------------------------------------------#

# The basic strategy for computing log(x) is to set r = log(x) and use
# Newton's method to solve the equation exp(r) = x. We set the initial
# value r_0 to math.log(x) and then iterate r_{n+1} = r_n + exp(-r_n) - 1
# until convergence. As with square roots, we increase the working
# precision dynamically during the process so that only one full-precision
# evaluation of exp is required.

# log(x) is small for most inputs, so the r values can safely be
# computed using fixed-point arithmetic. However, when x has a very
# large or small exponent, we can improve performance through the
# normalization log(t * 2**n) = log(t) + n*log(2), choosing n such
# that 0.5 <= t <= 1 (for example).

# There are some caveats: if x is extremely close to 1, the working
# precision must be increased to maintain high relative precision in the
# output (alternatively, the series approximation for log(1+x) could
# be used in that case).

# This function performs the Newton iteration using fixed-point
# arithmetic. x is assumed to have magnitude ~= 1
def log_newton(x, prec):
    extra = 10
    # 50-bit approximation
    fx = math.log(x) - 0.69314718055994529*prec
    r = int(fx * 2.0**50)
    prevp = 50
    for p in giant_steps(50, prec+extra):
        rb = lshift(r, p-prevp)

        # Parameters for exponential series
        r = int(2 * p**0.4)
        exp_extra = r + 10

        e = exp_series((-rb) << exp_extra, p + exp_extra, r)
        r = rb + ((rshift(x, prec-p)*e)>>(p + exp_extra)) - (1 << p)
        prevp = p
    return r >> extra

def flog(x, prec, rnd=round_fast):
    """Compute the natural logarithm of a positive raw mpf x."""
    sign, man, exp, bc = x
    if not man:
        if x == fzero:
            return fnan
        if x == finf:
            return finf
        return fnan
    if sign:
        raise ComplexResult("logarithm of a negative number")
    if x == fone:
        return fzero
    bc_plus_exp = bc + exp
    # Estimated precision needed for log(t) + n*log(2)
    prec2 = prec + int(math.log(1+abs(bc_plus_exp), 2)) + 15
    # Watch out for the case when x is very close to 1
    if -1 < bc_plus_exp < 2:
        near_one = fabs(fsub(x, fone, 53), 53)
        if near_one == 0:
            return fzero
        # estimate how close
        prec2 += -(near_one[2]) - bitcount(abs(near_one[1]))
    # Separate mantissa and exponent, calculate, join parts
    t = rshift(man, bc-prec2)
    l = log_newton(t, prec2)
    a = bc_plus_exp * log2_fixed(prec2)
    return from_man_exp(l+a, -prec2, prec, rnd)


#----------------------------------------------------------------------------#
#                                                                            #
#                          Trigonometric functions                           #
#                                                                            #
#----------------------------------------------------------------------------#

# We compute sin(x) around 0 from its Taylor series, and cos(x) around 0
# from sqrt(1-sin(x)**2). This way we can simultaneously compute sin and
# cos, which are often needed together (e.g. for the tangent function or
# the complex exponential), with little extra cost compared to computing
# just one of them. The main reason for computing sin first (and not sin
# from cos) is to obtain high relative accuracy for x extremely close to
# 0, where the operation sqrt(1-cos(x)**2) can cause huge cancellations.

# For any value of x, we can reduce it to the interval A = [-pi/4, pi/4]
# (where the Taylor series converges quickly) by translations, changing
# signs, and switching the roles of cos and sin:

#    A : sin(x) = sin(x)           cos(x) = cos(x)
#    B : sin(x) = cos(x-pi/2)      cos(x) = -sin(x-pi/2)
#    C : sin(x) = -sin(x-pi)       cos(x) = -cos(x-pi)
#    D : sin(x) = -cos(x-3*pi/2)   cos(x) = sin(x-3*pi/2)

# |     A      |      B     |      C     |     D     |
# v            v            v            v           v
# 
#    1 |  ____   ..........                            ____
#      |      _..          ..                        __
#      |      . __           .                     __
#      |    ..    _           ..                  _
#      |   .       __           .               __
# -----| -.----------_-----------.-------------_-----------
#      | .            _           ..          _           .
#      |               __           .       __           .
#      |                 _           ..    _           ..
#      |                  __           . __           .
#      |                    __         _..          ..
#   -1 |                      _________   ..........
#       0                       pi                     2*pi


def sin_taylor(x, prec):
    x2 = (x*x) >> prec
    s = a = x
    k = 3
    while a:
        a = ((a * x2) >> prec) // (k*(1-k))
        s += a
        k += 2
    return s

def cos_taylor(x, prec):
    x2 = (x*x) >> prec
    a = c = (1<<prec)
    k = 2
    while a:
        a = ((a * x2) >> prec) // (k*(1-k))
        c += a
        k += 2
    return c

# TODO: what we really want is the general version of Ziv's algorithm
def clamp1(x, prec, rnd):
    """Ensure that cos and sin values round away from -1 and 1 when
    using directed rnd."""
    if rnd is round_nearest or x[1] != 1:
        return x
    if rnd in (round_down, round_floor):
        if x == fone: return (0, (1<<prec)-1, -prec, prec)
    if rnd in (round_down, round_ceiling):
        if x == fnone: return (1, (1<<prec)-1, -prec, prec)
    return x

def cos_sin(x, prec, rnd=round_fast, nt=2):
    """Simultaneously compute (cos(x), sin(x)) for real x
    for nt=2; for nt=0 (1) compute only cos (sin)
    """

    sign, man, exp, bc = x

    if not man:
        if exp:
            c, s = (fnan, fnan)
        else:
            c, s = fone, fzero
        if   nt == 2:  return c, s
        elif nt == 0: return c
        else: return s

    magnitude = bc + exp
    abs_mag = abs(magnitude)

    # Very close to 0
    if magnitude < -prec:
        # Essentially exact
        if rnd is round_nearest:
            c, s = fone, fpos(x, prec, rnd)
            if   nt == 2: return c, s
            elif nt == 0: return c
            else: return s
        # Magic for interval arithmetic
        # cos(x) lies between 1 and 1-eps(1)/2
        if rnd in (round_up, round_ceiling):
            c = fone
        else:
            c = (0, (1<<prec)-1, -prec, prec)
        # sin(x) lies between x and x-eps(x)/2
        if rnd in (round_up, [round_ceiling, round_floor][sign]):
            s = fpos(x, prec, rnd)
        elif sign:
            s = fadd(x, (0, 1, magnitude-prec-4, 1), prec, rnd)
        else:
            s = fadd(x, (1, 1, magnitude-prec-4, 1), prec, rnd)
        if   nt == 2: return c, s
        elif nt == 0: return c
        else:         return s

    wp = wp1 = prec + 22

    while 1:
        # Reduce modulo pi/4
        rp = wp + abs_mag
        a = to_fixed(x, rp)
        pi = pi_fixed(rp)
        pi4 = pi >> 2
        pi2 = pi >> 1
        n, rx = divmod(a+pi4, pi2)
        rx -= pi4
        rx >>= abs_mag

        # If we're close to a root, we have to increase the
        # fixed-point precision to obtain full relative accuracy
        if abs(rx >> (wp1 - 8)) < 10:
            wp += wp1 - bitcount(abs(rx))
        else:
            break

    case = n % 4
    if nt == 2:
        one = 1 << wp
        s = sin_taylor(rx, wp)
        c = sqrt_fixed(one - ((s*s)>>wp), wp)
        if   case == 1: c, s = -s, c
        elif case == 2: s, c = -s, -c
        elif case == 3: c, s = s, -c
        c = clamp1(from_man_exp(c, -wp, prec, rnd), prec, rnd)
        s = clamp1(from_man_exp(s, -wp, prec, rnd), prec, rnd)
        return c, s
    elif nt == 0:
        if   case == 0: c = cos_taylor(rx, wp)
        elif case == 1: c = -sin_taylor(rx, wp)
        elif case == 2: c = -cos_taylor(rx, wp)
        elif case == 3: c = sin_taylor(rx, wp)
        return clamp1(from_man_exp(c, -wp, prec, rnd), prec, rnd)
    else:
        if   case == 0: s = sin_taylor(rx, wp)
        elif case == 1: s = cos_taylor(rx, wp)
        elif case == 2: s = -sin_taylor(rx, wp)
        elif case == 3: s = -cos_taylor(rx, wp)
        return clamp1(from_man_exp(s, -wp, prec, rnd), prec, rnd)

def fcos(x, prec, rnd=round_fast):
    return cos_sin(x, prec, rnd, 0)

def fsin(x, prec, rnd=round_fast):
    return cos_sin(x, prec, rnd, 1)

def ftan(x, prec, rnd=round_fast):
    c, s = cos_sin(x, prec+6)
    return fdiv(s, c, prec, rnd)


#----------------------------------------------------------------------
# Hyperbolic functions
#

def sinh_taylor(x, prec):
    x2 = (x*x) >> prec
    s = a = x
    k = 3
    while a:
        a = ((a * x2) >> prec) // (k*(k-1))
        s += a
        k += 2
    return s

def cosh_sinh(x, prec, rnd=round_fast):
    """Simultaneously compute (cosh(x), sinh(x)) for real x"""
    sign, man, exp, bc = x
    if (not man) and exp:
        if x == finf: return (finf, finf)
        if x == fninf: return (finf, fninf)
        return fnan, fnan

    if sign:
        man = -man

    high_bit = exp + bc
    prec2 = prec + 6

    if high_bit < -3:
        # Extremely close to 0, sinh(x) ~= x and cosh(x) ~= 1
        # TODO: support directed rounding
        if high_bit < -prec-2:
            return (fone, fpos(x, prec, rnd))

        # Avoid cancellation when computing sinh
        # TODO: might be faster to use sinh series directly
        prec2 += (-high_bit) + 4

    # In the general case, we use
    #    cosh(x) = (exp(x) + exp(-x))/2
    #    sinh(x) = (exp(x) - exp(-x))/2
    # and note that the exponential only needs to be computed once.
    ep = fexp(x, prec2)
    em = fdiv(fone, ep, prec2)
    ch = fshift(fadd(ep, em, prec, rnd), -1)
    sh = fshift(fsub(ep, em, prec, rnd), -1)
    return ch, sh

def fcosh(x, prec, rnd=round_fast):
    """Compute cosh(x) for a real argument x"""
    return cosh_sinh(x, prec, rnd)[0]

def fsinh(x, prec, rnd=round_fast):
    """Compute sinh(x) for a real argument x"""
    return cosh_sinh(x, prec, rnd)[1]

def ftanh(x, prec, rnd=round_fast):
    """Compute tanh(x) for a real argument x"""
    ch, sh = cosh_sinh(x, prec+6)
    return fdiv(sh, ch, prec, rnd)


#----------------------------------------------------------------------
# Inverse tangent
#

"""
Near x = 0, use atan(x) = x - x**3/3 + x**5/5 - ...
Near x = 1, use atan(x) = y/x * (1 + 2/3*y + 2*4/3/5*y**2 + ...)
where y = x**2/(1+x**2).

TODO: these series are not impressively fast. It is probably better
to calculate atan from tan, using Newton's method or even the
secant method.
"""

def atan_taylor(x, prec, rnd=round_fast):
    sign, man, exp, bc = x
    assert not sign
    # Increase absolute precision when extremely close to 0
    bc = bitcount(man)
    diff = -(bc + exp)
    prec2 = prec
    if diff > 10:
        if 3*diff - 4 > prec:  # x**3 term vanishes; atan(x) ~x
            return from_man_exp(man, exp, prec, rnd)
        prec2 = prec + diff
    prec2 += 15  # XXX: better estimate for number of guard bits
    x = to_fixed(x, prec2)
    x2 = (x*x)>>prec2; one = 1<<prec2; s=a=x
    for n in xrange(1, 1000000):
        a = (a*x2) >> prec2
        s += a // ((-1)**n * (n+n+1))
        if -100 < a < 100:
            break
    return from_man_exp(s, -prec2, prec, rnd)

def atan_euler(x, prec, rnd=round_fast):
    prec2 = prec + 15
    x = to_fixed(x, prec2)
    one = 1<<prec2; x2 = (x*x)>>prec2; y=(x2<<prec2)//(one+x2)
    s = a = one
    for n in xrange(1, 1000000):
        a = ((a*y)>>prec2) * (2*n) // (2*n+1)
        if a < 100:
            break
        s += a
    return from_man_exp(y*s//x, -prec2, prec, rnd)

_cutoff_1 = (0, 5, -3, 3)   # ~0.6
_cutoff_2 = (0, 3, -1, 2)   # 1.5

def fatan(x, prec, rnd=round_fast):
    sign, man, exp, bc = x
    if not man:
        if x == fzero: return fzero
        if x == finf: return fshift(fpi(prec), -1)
        if x == fninf: return fneg(fshift(fpi(prec), -1))
        return fnan
    if sign:
        return fneg(fatan(fneg(x), prec, rnd))
    if fcmp(x, _cutoff_1) < 0:
        return atan_taylor(x, prec, rnd)
    if fcmp(x, _cutoff_2) < 0:
        return atan_euler(x, prec, rnd)
    # For large x, use atan(x) = pi/2 - atan(1/x)
    if x[2] > 10*prec:
        pi = fpi(prec, rnd)
        pihalf = fshift(pi, -1)
    else:
        pi = fpi(prec+4)
        pihalf = fshift(pi, -1)
        t = fatan(fdiv(fone, x, prec+4), prec+4)
        return fsub(pihalf, t, prec, rnd)

def fatan2(y, x, prec, rnd=round_fast):
    xsign, xman, xexp, xbc = x
    ysign, yman, yexp, ybc = y
    if ysign:
        return fneg(fatan2(fneg(y), x, prec, rnd))
    if not xman:
        if not yman:
            return fzero
        return fshift(fpi(prec, rnd), -1)
    tquo = fatan(fdiv(y, x, prec+4), prec+4)
    if xsign:
        return fadd(fpi(prec+4), tquo, prec, rnd)
    else:
        return fpos(tquo, prec, rnd)

def fasin(x, prec, rnd=round_fast):
    # asin(x) = 2*atan(x/(1+sqrt(1-x**2)))
    sign, man, exp, bc = x
    if bc+exp > 0 and x not in (fone, fnone):
        raise ComplexResult("asin(x) is real only for -1 <= x <= 1")
    wp = prec + 15
    a = fmul(x, x)
    b = fadd(fone, fsqrt(fsub(fone, a, wp), wp), wp)
    c = fdiv(x, b, wp)
    return fshift(fatan(c, prec, rnd), 1)

def facos(x, prec, rnd=round_fast):
    # acos(x) = 2*atan(sqrt(1-x**2)/(1+x))
    sign, man, exp, bc = x
    if bc + exp > 0:
        if x not in (fone, fnone):
            raise ComplexResult("acos(x) is real only for -1 <= x <= 1")
        if x == fnone:
            return fpi(prec, rnd)
    wp = prec + 15
    a = fmul(x, x)
    b = fsqrt(fsub(fone, a, wp), wp)
    c = fdiv(b, fadd(fone, x, wp), wp)
    return fshift(fatan(c, prec, rnd), 1)

def fasinh(x, prec, rnd=round_fast):
    # asinh(x) = log(x+sqrt(x**2+1))
    wp = prec + 15
    q = fsqrt(fadd(fmul(x,x), fone, wp), wp)
    return flog(fadd(x, q, wp), prec, rnd)

def facosh(x, prec, rnd=round_fast):
    # acosh(x) = log(x+sqrt(x**2-1))
    wp = prec + 15
    if fcmp(x, fone) == -1:
        raise ComplexResult("acosh(x) is real only for x >= 1")
    q = fsqrt(fadd(fmul(x,x), fnone, wp), wp)
    return flog(fadd(x, q, wp), prec, rnd)

def fatanh(x, prec, rnd=round_fast):
    # atanh(x) = log((1+x)/(1-x))/2
    sign, man, exp, bc = x
    mag = bc + exp
    if mag > 0:
        raise ComplexResult("atanh(x) is real only for -1 < x < 1")
    wp = prec + 15
    a = fadd(x, fone, wp)
    b = fsub(fone, x, wp)
    return fshift(flog(fdiv(a, b, wp), prec, rnd), -1)

