from __future__ import print_function, division

import math

from collections import (
    Counter,
    OrderedDict,
)


def count_digits(n, base=10):
    """
    Returs an ordered dict. frequency counter for the digits of a decimal
    integer ``n`` representing an integer in a base ``base``, where the
    base digits are mapped to the positive decimal integers, e.g. the
    base ``11`` integer ``10321``, with digits ``10``, ``3``, ``2``, ``1``
    (from most significant to least), is equal to the decimal integer
    ``13696``, as ``10 * 11**3 + 3 * 11**2 + 2 * 11**1 + 1 * 11**0 = 13696``.

    For binary, octal or hexadecimal bases the integer ``n``
    can also be specified as a numeric literal in those bases, e.g. ``0xF321``.

    Examples
    ========

    >>> sympy.ntheory import count_digits

    >>> count_digits(1903413094, 10)
    OrderedDict([(0, 2), (1, 2), (3, 2), (4, 2), (9, 2)])

    >>> count_digits(122333, 10)
    OrderedDict([(1, 1), (2, 2), (3, 3)])

    >>> count_digits(0b1001001), 2)
    OrderedDict([(0, 4), (1, 3)])

    >>> count_digits(0xF321, 16)
    OrderedDict([(1, 1), (2, 1), (3, 1), (15, 1)])

    The base 11 number 10321, with digits 10, 3, 2, 1 (from most significant digit)
    >>> count_digits(13696, 11)
    OrderedDict([(1, 1), (2, 1), (3, 1), (10, 1)]

    The base 100 number 990512, with digits 99, 0, 51, 2 (from most significant digit)
    >>> count_digits(99005102, 100)
    OrderedDict([(0, 1), (2, 1), (51, 1), (99, 1)])

    """
    m = math.floor(math.log(n, base)) + 1

    return OrderedDict(
        sorted(
            Counter(((n // base ** k) % base for k in range(m))).items(),
            key=lambda t: t[0]
        )
    )


def is_palindromic_integer(n, base=10):
    """
    Checks whether a given positive integer ``n``, specified either as a
    numeric literal (binary, octal, decimal or hexadecimal) or a string, has
    a symmetrical digit sequence, i.e. the sequence of its digits is the same
    left to right vs right to left.

    Examples
    ========

    >>> is_palindromic_integer(1, 10)
    True

    >>> is_palindromic_integer(33, 10)
    True

    >>> is_palindromic_integer(919, 10)
    True

    >>> is_palindromic_integer(0b101, 2)
    True

    >>> is_palindromic_integer(0b01010010,2)
    False

    >>> is_palindromic_integer(0o7610167, 8)
    True

    >>> is_palindromic_integer(0o7611167, 8)
    False

    >>> is_palindromic_integer(0xFA0100AF, 16)
    False

    >>> is_palindromic_integer(0xFA0110AF, 16)
    True

    """
    if base < 2:
        raise ValueError('The base must be at least 2')
    s = str(n)
    m = len(s)

    return s[:int(m / 2)] == (s[int(m / 2):] if not m % 2 else s[int(m / 2) + 1:])[::-1]
