""" This module is an alternative implementation of the gmpy2 function.

In sympy, gmpy2 is optional. This means that we need to provide
some alternative function for environments where gmpy2 is not installed.
On the other hand, for developers, we want to provide functions without making
them aware of the differences. `sympy.external.gmpy` is provided
with such an intention. That is, when gmpy2 is installed,
the functions of gmpy2 are used, and when gmpy2 is not installed,
the functions defined in this module are used. The choice of which
function to use is handled by `sympy.external.gmpy`,
which implements the alternative functions that may be selected.
Therefore, it is not expected that functions in this module
will be called directly by the user.
"""

import math
import mpmath.libmp as mlib


def factorial_python(x):
    return int(mlib.ifac(x))


def sqrt_python(x):
    return int(mlib.isqrt(x))


def is_square_python(x):
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
    return mlib.sqrtrem(x)[1] == 0


def sqrtrem_python(x):
    return tuple(int(r) for r in mlib.sqrtrem(x))


def invert_python(x, m):
    """ Return y such that x*y == 1 modulo m.

    Uses ``math.pow`` but reproduces the behaviour of ``gmpy2.invert``
    which raises ZeroDivisionError if no inverse exists.
    """
    try:
        return pow(x, -1, m)
    except ValueError:
        raise ZeroDivisionError("invert() no inverse exists")


def legendre_python(x, y):
    """ Return Legendre symbol (x / y).

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


def jacobi_python(x, y):
    """ Return Jacobi symbol (x / y)."""
    if y <= 0 or not y % 2:
        raise ValueError("y should be an odd positive integer")
    x %= y
    if not x:
        return int(y == 1)
    if y == 1 or x == 1:
        return 1
    if math.gcd(x, y) != 1:
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
