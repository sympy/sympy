from __future__ import print_function, division

from sympy import Rational
from fractions import gcd
from sympy.core.numbers import igcdex
from math import log

def egyptian_fraction(r, algorithm="Greedy"):
    """
    Return the list of denominators of an Egyptian fraction
    expansion [1] of the said rational `r`.

    Parameters
    ==========

    r : Rational
        a rational number
    algorithm : { "Greedy", "Graham Jewett", "Takenouchi", "Golomb"}, optional
        Denotes the algorithm to be used (the default is "Greedy").

    Examples
    ========

    >>> from sympy import Rational
    >>> from sympy.ntheory.egyptian_fraction import egyptian_fraction
    >>> egyptian_fraction(Rational(3, 7))
    [3, 11, 231]
    >>> egyptian_fraction(Rational(3, 7), "Graham Jewett")
    [7, 8, 9, 56, 57, 72, 3192]
    >>> egyptian_fraction(Rational(3, 7), "Takenouchi")
    [4, 7, 28]
    >>> egyptian_fraction(Rational(3, 7), "Golomb")
    [3, 15, 35]

    See Also
    ========

    sympy.core.numbers.Rational

    Notes
    =====

    Currently the following algorithms are supported:

    1) Greedy Algorithm
        Also called the Fibonacci-Sylvester algorithm [2].
        At each step, extract the largest unit fraction less
        than the target and replace the target with the remainder.
        Limits
        ======
        a) Given *p/q* in lowest terms, generates an expansion of maximum
            length *p*. Even as the numerators get large, the number of
            terms is seldom more than a handful.

        b) Uses minimal memory.

        c) Only potential problem is that the terms can blow up (standard
            examples of this are 5/121 and 31/311); however, this isn't really
            a problem in Python until the numbers reach obscenely large numbers
            of digits.  The denominator is at most squared at each step
            (doubly-exponential growth) and typically exhibits singly-exponential
            growth.

    2) Graham Jewett Algorithm
        The algorithm suggested by the result of Graham and Jewett.
        Note that this has a tendency to blow up: the length of the
        resulting expansion is always 2**(x/gcd(x, y)) - 1 .  See [3].

    3) Takenouchi Algorithm
        The algorithm suggested by Takenouchi (1921).
        Differs from the Graham-Jewett algorithm only in the handling
        of duplicates.  See [3].

    4) Golomb's Algorithm
        Developed by Golomb, this algorithm is significantly faster and has the
        advantage of not failing when x > y. Produces x or fewer terms, none
        greater than y*(y-1).

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Egyptian_fraction
    .. [2] https://en.wikipedia.org/wiki/Greedy_algorithm_for_Egyptian_fractions
    .. [3] http://www.ics.uci.edu/~eppstein/numth/egypt/conflict.html

    """
    x, y = r.as_numer_denom()

    if algorithm == "Greedy":
        return egypt_greedy(x, y)
    elif algorithm == "Graham Jewett":
        return egypt_graham_jewett(x, y)
    elif algorithm == "Takenouchi":
        return egypt_takenouchi(x, y)
    elif algorithm == "Golomb":
        return egypt_golomb(x, y)
    else:
        raise ValueError("Entered invalid algorithm")


def egypt_greedy(x, y):
    if x == 1:
        return [y]
    else:
        a = (-y) % (x)
        b = y*(y//x + 1)
        c = gcd(a, b)
        if c > 1:
            num, denom = a//c, b//c
        else:
            num, denom = a, b
        return [y//x + 1] + egypt_greedy(num, denom)


def egypt_graham_jewett(x, y):
    l = [y] * x

    # l is now a list of integers whose reciprocals sum to x/y.
    # we shall now proceed to manipulate the elements of l without
    # changing the reciprocated sum until all elements are unique.

    while len(l) != len(set(l)):
        l.sort()  # so the list has duplicates. find a smallest pair
        for i in range(len(l) - 1):
            if l[i] == l[i + 1]:
                break
        # we have now identified a pair of identical
        # elements: l[i] and l[i + 1].
        # now comes the application of the result of graham and jewett:
        l[i + 1] = l[i] + 1
        # and we just iterate that until the list has no duplicates.
        l.append(l[i]*(l[i] + 1))
    return sorted(l)


def egypt_takenouchi(x, y):
    l = [y] * x
    while len(l) != len(set(l)):
        l.sort()
        for i in range(len(l) - 1):
            if l[i] == l[i + 1]:
                break
        k = l[i]
        if k % 2 == 0:
            l[i] = l[i] // 2
            del l[i + 1]
        else:
            l[i], l[i + 1] = (k + 1)//2, k*(k + 1)//2
    return sorted(l)

def modinv(a, m):
    # Returns the multiplicative inverse of a modulo m
    x, y, g = igcdex(a, m)
    if g != 1:
        raise Exception("%d has no invesrse modulo %d" % a, m)
    else:
        return x % m

def egypt_golomb(x, y):
    expansion = []
    div = gcd(x, y)
    x, y = x/div, y/div
    x2, y2 = x, y
    while True:
        if x2 == 1:
            expansion.append(int(y2))
            break
        x3 = modinv(x2, y2)
        z = (x2 * x3 - 1) / y2
        expansion.append(int(x3 * y2))
        y2, x2 = x3, z
    return sorted(expansion)
