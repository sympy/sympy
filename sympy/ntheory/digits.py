from __future__ import print_function, division

from collections import defaultdict

from sympy.core.compatibility import as_int
from sympy.utilities.iterables import multiset, is_palindromic as _palindromic


def digits(n, b=10):
    """
    Return a list of the digits of n in base b. The first element in the list
    is b (or -b if n is negative).

    Examples
    ========

    >>> from sympy.ntheory.digits import digits
    >>> digits(35)
    [10, 3, 5]
    >>> digits(27, 2)
    [2, 1, 1, 0, 1, 1]
    >>> digits(65536, 256)
    [256, 1, 0, 0]
    >>> digits(-3958, 27)
    [-27, 5, 11, 16]
    """

    b = as_int(b)
    n = as_int(n)
    if b <= 1:
        raise ValueError("b must be >= 2")
    else:
        x, y = abs(n), []
        while x >= b:
            x, r = divmod(x, b)
            y.append(r)
        y.append(x)
        y.append(-b if n < 0 else b)
        y.reverse()
        return y


def count_digits(n, base=10):
    """
    Return a dictionary whose keys are the digits of ``n`` in the
    given base, with keys indicating the digits appearing in the
    number and values indicating how many times that digit appeared.

    Examples
    ========

    >>> from sympy.ntheory import count_digits, digits

    >>> count_digits(1111339)
    {1: 4, 3: 2, 9: 1}

    The digits returned are always represented in base-10
    but the number itself can be entered in any format that is
    understood by Python; the base of the number can also be
    given if it is different than 10:

    >>> n = 0xFA; n
    250
    >>> count_digits(_)
    {0: 1, 2: 1, 5: 1}
    >>> count_digits(n, 16)
    {10: 1, 15: 1}

    The default dictionary will return a 0 for any digit that did
    not appear in the number. For example, which digits appear 7
    times in ``77!``:

    >>> from sympy import factorial
    >>> c77 = count_digits(factorial(77))
    >>> [i for i in range(10) if c77[i] == 7]
    [1, 3, 7, 9]
    """
    return defaultdict(int, multiset(digits(n, base)[1:]).items())


def is_palindromic(n, base=10):
    """return True if ``n`` is the same when read from left to right
    or right to left in the given base.

    Examples
    ========

    >>> from sympy.ntheory import is_palindromic

    >>> all(is_palindromic(i) for i in (-11, 1, 22, 121))
    True

    The ``base`` argument allows you to test numbers in other
    bases. For example, 88 is palindromic in base-10 but not
    in base-8:

    >>> is_palindromic(88, 8)
    False

    On the other hand, a number can be palindromic in base-8 but
    not in base-10:

    >>> is_palindromic(0o121)
    False
    >>> 0o121
    81

    Or it might be palindromic in both bases:

    >>> is_palindromic(121, 8) and is_palindromic(121)
    True
    >>> 0o171 == 121
    True

    """
    try:
        n = as_int(n)
    except ValueError:
        n = int(n)
    base = as_int(base)

    if base < 2:
        raise ValueError('The base must be at least 2')

    if n < 0:
        n = -n

    _str = {2: bin, 8: oct, 16: hex}.get(base, '')
    if _str:
        s = _str(n)[2:]  # ignore sign and 2 char base indicator
    else:
        s = str(n)
    return _palindromic(s)
