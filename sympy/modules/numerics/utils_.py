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

def isqrt(y):
    """Calculate floor square root of an integer"""
    # Start from regular floating-point estimate. Rewrite as
    # sqrt(y)=2**(log_2(y)/2-n)*2**n to avoid overflow.
    if y == 0:
        return 0
    lg = _clog(y, 2)/2
    n = max(int(lg)-52, 0)
    guess = int(2.0**(lg-n)) << n
    # Newton iteration
    xprev, x = -1, guess
    while abs(x - xprev) > 1:
        xprev, x = x, (x + y//x)>>1
    while x*x > y:
        x -= 1
    return x

def make_fixed(x, prec):
    """Convert a Float to a fixed-point big integer"""
    offset = x.exp+prec
    if offset >= 0:
        return x.man << offset
    else:
        return x.man >> (-offset)
