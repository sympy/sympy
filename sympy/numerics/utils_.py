from math import log as _clog

_bc_table = [0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4]

def bitcount(n):
    """Give position of the highest set bit in an integer"""
    if n == 0: return 0
    if n < 0: n = -n
    # math.log gives a good estimate, and never overflows, but
    # is not always exact. Subtract 2 to underestimate, then
    # count remaining bits by table lookup
    bc = max(0, int(_clog(n, 2)) - 2)
    return bc + _bc_table[n >> bc]

def trailing_zeros(n):
    """Count trailing zero bits in an integer."""
    if n & 1: return 0
    if n == 0: return 0
    if n < 0: n = -n
    t = 0
    while not n & 0xffffffffffffffff: n >>= 64; t += 64
    while not n & 0xff: n >>= 8; t += 8
    while not n & 1: n >>= 1; t += 1
    return t

def make_fixed(x, prec):
    """Convert a Float to a fixed-point big integer"""
    offset = x.exp+prec
    if offset >= 0:
        return x.man << offset
    else:
        return x.man >> (-offset)

def bin_to_radix(x, xbits, base, bdigits):
    return x * (base**bdigits) >> xbits

_numerals = '0123456789abcdefghijklmnopqrstuvwxyz'

def small_numeral(n, base=10):
    # Calculate numeral of n*(base**digits) in the given base
    if base == 10:
        return str(n)
    digs = []
    while n:
        n, digit = divmod(n, base)
        digs.append(_numerals[digit])
    return "".join(digs[::-1])

# TODO: speed up for bases 2, 4, 8, 16, ...
def fixed_to_str(x, base, digits):
    if digits < 789:
        return small_numeral(x, base)
    half = (digits // 2) + (digits & 1)
    if "verbose" in global_options and half > 50000:
        print "  dividing..."
    A, B = divmod(x, base**half)
    ad = fixed_to_str(A, base, half)
    bd = fixed_to_str(B, base, half).rjust(half, "0")
    return ad + bd

global_options = {}
