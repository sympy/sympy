# sympy.external.ntheory
#
# This module provides pure Python implementations of some number theory
# functions that are alternately used from gmpy2 if it is installed.

import sys
import math

import mpmath.libmp as mlib


def factorial(x):
    """Return x!."""
    return int(mlib.ifac(int(x)))


def sqrt(x):
    """Integer square root of x."""
    return int(mlib.isqrt(int(x)))


def sqrtrem(x):
    """Integer square root of x and remainder."""
    s, r = mlib.sqrtrem(int(x))
    return (int(s), int(r))


if sys.version_info[:2] >= (3, 9):
    # As of Python 3.9 these can take multiple arguments
    gcd = math.gcd
    lcm = math.lcm

else:
    # Until python 3.8 is no longer supported
    from functools import reduce


    def gcd(*args):
        """gcd of multiple integers."""
        return reduce(math.gcd, args, 0)


    def lcm(*args):
        """lcm of multiple integers."""
        if 0 in args:
            return 0
        return reduce(lambda x, y: x*y//math.gcd(x, y), args, 1)


def is_square(x):
    """Return True if x is a square number."""
    if x < 0:
        return False

    # Note that the possible values of y**2 % n for a given n are limited.
    # For example, when n=4, y**2 % n can only take 0 or 1.
    # In other words, if x % 4 is 2 or 3, then x is not a square number.
    # Mathematically, it determines if it belongs to the set {y**2 % n},
    # but implementationally, it can be realized as a logical conjunction
    # with an n-bit integer.
    # see https://mersenneforum.org/showpost.php?p=110896
    # def magic(n):
    #     s = {y**2 % n for y in range(n)}
    #     s = set(range(n)) - s
    #     return sum(1 << bit for bit in s)
    # >>> print(hex(magic(128)))
    # 0xfdfdfdedfdfdfdecfdfdfdedfdfcfdec
    # >>> print(hex(magic(99)))
    # 0x5f6f9ffb6fb7ddfcb75befdec
    # >>> print(hex(magic(91)))
    # 0x6fd1bfcfed5f3679d3ebdec
    # >>> print(hex(magic(85)))
    # 0xdef9ae771ffe3b9d67dec
    if 0xfdfdfdedfdfdfdecfdfdfdedfdfcfdec & (1 << (x & 127)):
        return False  # e.g. 2, 3
    m = x % 765765 # 765765 = 99 * 91 * 85
    if 0x5f6f9ffb6fb7ddfcb75befdec & (1 << (m % 99)):
        return False  # e.g. 17, 68
    if 0x6fd1bfcfed5f3679d3ebdec & (1 << (m % 91)):
        return False  # e.g. 97, 388
    if 0xdef9ae771ffe3b9d67dec & (1 << (m % 85)):
        return False  # e.g. 793, 1408
    return mlib.sqrtrem(int(x))[1] == 0


def invert(x, m):
    """Modular inverse of x modulo m.

    Returns y such that x*y == 1 mod m.

    Uses ``math.pow`` but reproduces the behaviour of ``gmpy2.invert``
    which raises ZeroDivisionError if no inverse exists.
    """
    try:
        return pow(x, -1, m)
    except ValueError:
        raise ZeroDivisionError("invert() no inverse exists")


def legendre(x, y):
    """Legendre symbol (x / y).

    Following the implementation of gmpy2,
    the error is raised only when y is an even number.
    """
    if y <= 0 or not y % 2:
        raise ValueError("y should be an odd prime")
    x %= y
    if not x:
        return 0
    if pow(x, (y - 1) // 2, y) == 1:
        return 1
    return -1


def jacobi(x, y):
    """Jacobi symbol (x / y)."""
    if y <= 0 or not y % 2:
        raise ValueError("y should be an odd positive integer")
    x %= y
    if not x:
        return int(y == 1)
    if y == 1 or x == 1:
        return 1
    if gcd(x, y) != 1:
        return 0
    j = 1
    while x != 0:
        while x % 2 == 0 and x > 0:
            x >>= 1
            if y % 8 in [3, 5]:
                j = -j
        x, y = y, x
        if x % 4 == y % 4 == 3:
            j = -j
        x %= y
    return j


def kronecker(x, y):
    """Kronecker symbol (x / y)."""
    if gcd(x, y) != 1:
        return 0
    if y == 0:
        return 1
    sign = -1 if y < 0 and x < 0 else 1
    y = abs(y)
    # We want to calculate s = trailing(y)
    s = 0
    while y % 2 == 0:
        y >>= 1
        s += 1
    if s % 2 and x % 8 in [3, 5]:
        sign = -sign
    return sign * jacobi(x, y)


def iroot(y, n):
    if y < 0:
        raise ValueError("y must be nonnegative")
    if n < 1:
        raise ValueError("n must be positive")
    if y in (0, 1):
        return y, True
    if n == 1:
        return y, True
    if n == 2:
        x, rem = mlib.sqrtrem(y)
        return int(x), not rem
    if n >= y.bit_length():
        return 1, False
    # Get initial estimate for Newton's method. Care must be taken to
    # avoid overflow
    try:
        guess = int(y**(1./n) + 0.5)
    except OverflowError:
        exp = math.log(y, 2)/n
        if exp > 53:
            shift = int(exp - 53)
            guess = int(2.0**(exp - shift) + 1) << shift
        else:
            guess = int(2.0**exp)
    if guess > 2**50:
        # Newton iteration
        xprev, x = -1, guess
        while 1:
            t = x**(n - 1)
            xprev, x = x, ((n - 1)*x + y//t)//n
            if abs(x - xprev) < 2:
                break
    else:
        x = guess
    # Compensate
    t = x**n
    while t < y:
        x += 1
        t = x**n
    while t > y:
        x -= 1
        t = x**n
    return x, t == y
