from __future__ import print_function, division

import re

from collections import Counter


def is_pandigital(n, zeroless=False, freq='1+'):
    """
    Checks whether a positive integer ``n`` is pandigital with respect
    to the (decimal) base and has a digit frequency given by the
    string ``freq``. The optional ``zeroless`` argument can be used
    to indicate whether ``0`` should or should not be included in
    checking for the pandigital property, i.e. whether ``0`` should be
    included in the base digit set.

    Note: ``freq`` is an integer string with an optional ``+``
    suffix to indicate min. digit frequency - otherwise the digit
    frequency is taken to be exact.

    Some examples are given below.
    ::
        915286437, False, '1+'           ->   False
        915286437, True, '1+'            ->   True
        9152860437, False, '1+'          ->   True
        112233445566778899, False, '2'   -> False
        112233445566778899, True, '2'    -> True
        11223344556677889900, False, '2' -> True
    """
    digits = Counter(str(n))

    _dfreq, dfreq_min = re.match(r'(\d+)(\+)?', freq).groups()
    _dfreq = int(_dfreq)

    base = set(Counter(str(1234567890)).keys()).difference(['0'] if zeroless else [])

    return base.issubset(digits.keys()) and (
        min(digits.values()) == max(digits.values()) == _dfreq if dfreq_min != '+' else
        min(digits.values()) >= _dfreq
    )


def is_palindromic_integer(n):
    """
    Checks whether a given positive integer ``n`` has the property that the
    sequence of its digits is the same left to right as it is right to left.
    Examples below.
    ::

        1 -> True
        33 -> True
        919 - > True
        1234321 -? True
        14987778941 -> True
        14987678941 -> False
    """
    s = str(n)
    m = len(s)

    return s[:int(m / 2)] == (s[int(m / 2):] if not m % 2 else s[int(m / 2) + 1:])[::-1]
