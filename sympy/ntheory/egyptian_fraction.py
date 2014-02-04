from __future__ import print_function, division

from sympy import Rational
from fractions import gcd


def egyptian_fraction(rat, algorithm="Greedy"):
    """
    Return the denominators of an Egtyptian fraction expansion of
    the said rational.

    The egyptian_fraction function takes rat and algorithm as input.
    'rat' is of the form Rational(n, d) where n is numerator and
    d is denominator. 'algorithm' denotes the algorithm to be
    used. By default, the 'algorithm' field can be left empty and
    the greedy algorithm will be used. 'algorithm' supports
    values : "Greedy", "Graham Jewett", and "Takenouchi".

    Currently the following algorithms are supported :

    1) Greedy Algorithm :
    Greedy algorithm for Egyptian fraction expansion
    also called the Fibonacci-Sylvester algorithm.
    At each step, extract the largest unit fraction less
    than the target and replace the target with the remainder.

    2) Graham Jewett Algorithm :
    The algorithm suggested by the result of Graham and Jewett.
    Note that this has a tendency to blow up: the length of the
    resulting expansion is always 2**(x/gcd(x, y)) - 1 .

    3) Takenouchi Algorithm :
    The algorithm suggested by Takenouchi (1921).
    Differs from the Graham-Jewett algorithm only in the handling
    of duplicates.

    Examples
    ========

    >>> from sympy import Rational
    >>> from sympy.ntheory.egyptian_fraction import egyptian_fraction
    >>> egyptian_fraction(Rational(3,7))
    [3, 11, 231]
    >>> egyptian_fraction(Rational(3, 7), "Greedy")
    [3, 11, 231]
    >>> egyptian_fraction(Rational(3, 7), "Graham Jewett")
    [7, 8, 9, 56, 57, 72, 3192]
    >>> egyptian_fraction(Rational(3, 7), "Takenouchi")
    [4, 7, 28]

    References
    =========

    .. [1] https://en.wikipedia.org/wiki/Greedy_algorithm_for_Egyptian_fractions
    .. [2] http://en.wikipedia.org/wiki/Egyptian_fraction#Modern_number_theory
    .. [3] http://www.ics.uci.edu/~eppstein/numth/egypt/conflict.html

    """
    x, y = rat.as_numer_denom()

    if algorithm == "Greedy":
        return egypt_greedy(x, y)
    elif algorithm == "Graham Jewett":
        return egypt_graham_jewett(x, y)
    elif algorithm == "Takenouchi":
        return egypt_takenouchi(x, y)
    else:
        raise ValueError("Entered Invalid Algorithm")


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
