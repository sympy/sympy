"""
Contents of this module:

* Integer and bit-level operations
* Miscellaneous utilities

"""

import math
import decimal

# Same as standard Python float
STANDARD_PREC = 53

LOG2_10 = math.log(10,2)  # 3.3219...


# All supported rounding modes. We define them as integer constants for easy
# management, but change __repr__ to give more information on inspection

class Rounding(int):

    def __new__(cls, level, name):
        a = int.__new__(cls, level)
        a.name = name
        return a

    def __repr__(self):
        return self.name

    __str__ = __repr__

ROUND_DOWN    = Rounding(1, 'ROUND_DOWN')
ROUND_UP      = Rounding(2, 'ROUND_UP')
ROUND_FLOOR   = Rounding(3, 'ROUND_FLOOR')
ROUND_CEILING = Rounding(4, 'ROUND_CEILING')
ROUND_HALF_UP = Rounding(5, 'ROUND_HALF_UP')
ROUND_HALF_DOWN = Rounding(6, 'ROUND_HALF_DOWN')
ROUND_HALF_EVEN = Rounding(7, 'ROUND_HALF_EVEN')


def giant_steps(start, target):
    """Generate a list of precisions ranging from 'start' to 'target'
    that doubles with each step. This is used by quadratically
    convergent iterations (that is, Newton iterations), where we want
    to keep the precision at the same level as the accuracy in each
    step to minimize work.

    For example, to find a sequence of precisions to reach 1000 bits
    starting from a 53-bit estimate, giant_steps(53, 1000) gives

        [64, 126, 251, 501, 1000]

    So, if we start Newton's method with a 53-bit accurate initial
    guess, the first iteration should be carried out at 64-bit
    precision, the second at 126-bit precision, and so on.

    Note the conservative rounding (1000 to 501, etc); this is used
    guard against unit errors in the last place."""
    L = [target]
    while L[-1] > start*2:
        L = L + [L[-1]//2 + 1]
    return L[::-1]


def rshift_quick(x, n):
    """For an integer x, calculate x >> n with the fastest (floor)
    rounding. Unlike the plain Python expression (x >> n), n is
    allowed to be negative, in which case a left shift is performed."""
    if n >= 0: return x >> n
    else:      return x << (-n)


def lshift_quick(x, n):
    """For an integer x, calculate x << n. Unlike the plain Python
    expression (x << n), n is allowed to be negative, in which case a
    right shift with default (floor) rounding is performed."""
    if n >= 0: return x << n
    else:      return x >> (-n)


def make_fixed(s, prec):
    """Convert a floating-point number to a fixed-point big integer"""
    man, exp, bc = s
    offset = exp + prec
    if offset >= 0:
        return man << offset
    else:
        return man >> (-offset)


def bitcount(n, log=math.log, table=(0,1,2,2,3,3,3,3,4,4,4,4,4,4,4,4)):
    """Give size of n in bits; i.e. the position of the highest set bit
    in n. If n is negative, the absolute value is used. The bitcount of
    zero is taken to be 0."""

    if not n: return 0
    if n < 0: n = -n

    # math.log gives a good estimate, and never overflows, but
    # is not always exact. Subtract 2 to underestimate, then
    # count remaining bits by table lookup
    bc = int(log(n, 2)) - 2
    if bc < 0:
        bc = 0
    return bc + table[n >> bc]


# from decimal.py -- faster for small precs
def bitcount2(n, correction = {
        '0': 4, '1': 3, '2': 2, '3': 2,
        '4': 1, '5': 1, '6': 1, '7': 1,
        '8': 0, '9': 0, 'a': 0, 'b': 0,
        'c': 0, 'd': 0, 'e': 0, 'f': 0}):
    if n < 0:
        n = -n
    hex_n = "%x" % n
    return 4*len(hex_n) - correction[hex_n[0]]


def trailing_zeros(n):
    """Count trailing zero bits in an integer. If n is negative, it is
    replaced by its absolute value."""
    if n & 1: return 0
    if not n: return 0
    if n < 0: n = -n
    t = 0
    while not n & 0xffffffffffffffff: n >>= 64; t += 64
    while not n & 0xff: n >>= 8; t += 8
    while not n & 1: n >>= 1; t += 1
    return t

#----------------------------------------------------------------------------
# Integer shifting with directed rounding
#

def round_floor(x, n):
    if not n or not x: return x
    if n < 0: return x << -n
    return x >> n

def round_ceiling(x, n):
    if not n or not x: return x
    if n < 0: return x << -n
    return -((-x) >> n)

def round_down(x, n):
    if not n or not x: return x
    if n < 0: return x << -n
    if x > 0:
        return x >> n
    else:
        return -((-x) >> n)

def round_up(x, n):
    if not n or not x: return x
    if n < 0: return x << -n
    if x > 0:
        return -((-x) >> n)
    else:
        return x >> n

def round_half_up(x, n):
    if not n or not x: return x
    if n < 0: return x << -n
    positive = x > 0
    if positive: t = x >> (n-1)
    else:        t = (-x) >> (n-1)
    if t & 1:
        if positive: return (t>>1)+1
        else:        return -((t>>1)+1)
    if positive: return t>>1
    else:        return -(t>>1)

def round_half_down(x, n):
    if not n or not x: return x
    if n < 0: return x << -n
    positive = x > 0
    if positive: t = x >> (n-1)
    else:        t = (-x) >> (n-1)
    if t & 1 and x & ((1<<(n-1))-1):
        if positive: return (t>>1)+1
        else:        return -((t>>1)+1)
    if positive: return t>>1
    else:        return -(t>>1)

def round_half_even(x, n):
    if not n or not x: return x
    if n < 0: return x << -n
    positive = x > 0
    if positive: t = x >> (n-1)
    else:        t = (-x) >> (n-1)
    if t & 1 and (t&2 or x & ((1<<(n-1))-1)):
        if positive: return (t>>1)+1
        else:        return -((t>>1)+1)
    if positive:  return t>>1
    else:         return -(t>>1)


rounding_functions = {
  ROUND_FLOOR : round_floor,
  ROUND_CEILING : round_ceiling,
  ROUND_DOWN : round_down,
  ROUND_UP : round_up,
  ROUND_HALF_UP : round_half_up,
  ROUND_HALF_DOWN : round_half_down,
  ROUND_HALF_EVEN : round_half_even
}


def rshift(x, n, rounding):
    """Shift x (a plain Python integer) n bits to the right (i.e.,
    calculate x/(2**n)), and round to the nearest integer in accordance
    with the specified rounding mode. The exponent n may be negative,
    in which case x is shifted to the left (and no rounding is
    necessary)."""
    return rounding_functions[rounding](x, n)



# TODO: speed up for bases 2, 4, 8, 16, ...

def bin_to_radix(x, xbits, base, bdigits):
    """
    Radix conversion for fixed-point numbers. That is, convert
    x * 2**xbits to floor(x * 10**bdigits).
    """
    return x * (base**bdigits) >> xbits

stddigits = '0123456789abcdefghijklmnopqrstuvwxyz'

def small_numeral(n, base=10, digits=stddigits):
    """
    Return the string numeral of a positive integer in an arbitrary
    base. Most efficient for small input.
    """
    if base == 10:
        return str(n)
    digs = []
    while n:
        n, digit = divmod(n, base)
        digs.append(digits[digit])
    return "".join(digs[::-1])

def numeral(n, base=10, size=0, digits=stddigits):
    """
    Represent the integer n as a string of digits in the given base.
    Recursive division is used to make this function about 3x faster
    than Python's str() for converting integers to decimal strings.

    The 'size' parameters specifies the number of digits in n; this
    number is only used to determine splitting points and need not
    be exact.
    """

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
