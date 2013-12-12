from __future__ import print_function, division

from sympy import Rational
from fractions import gcd
from random import randint

def egypt_greedy(Rational):
    """
    Greedy algorithm for Egyptian fraction expansion
    Also called the Fibonacci-Sylvester algorithm
    At each step, extract the largest unit fraction less
    than the target and replace the target with the remainder
    Inputs: x/y is the target fraction, not necessarily in lowest terms
    y is any natural number
    x is any integer strictly between 0 and y

    References
    =========

    - https://en.wikipedia.org/wiki/Greedy_algorithm_for_Egyptian_fractions

    Examples
    ========

    >>> from sympy import Rational
    >>> from sympy.ntheory.egyptian_fraction import egypt_greedy
    >>> egypt_greedy(Rational(7,12))
    [2, 12]
    >>> egypt_greedy(Rational(2,3))
    [2, 6]

    See Also
    =======

    egypt_graham_jewett(Rational) : Uses Graham and Jewett Algorithm
    egypt_takenouchi(Rational) : Uses Takenouchi Algorithm
    """

    x, y = int(Rational.p), int(Rational.q)
    if (x == 1):
        return [y]
    else:
        a = (-y) % (x)
        b = y * (int(y/x) + 1)
        c = gcd(a, b)
        if c > 1:
            num, denom = int(a/c), int(b/c)
        else:
            num, denom = a, b
        return [int(y/x) + 1] + egypt_greedy(Rational(num,denom))

def egypt_graham_jewett(Rational):
    """
    The algorithm suggested by the result of Graham and Jewett.
    Note that this has a tendency to blow up: the length of the resulting expansion
    is always 2**(x/gcd(x,y)) - 1,
    Same arguments as egypt_greedy(Rational).

    References
    =========

    - http://en.wikipedia.org/wiki/Egyptian_fraction#Modern_number_theory

    Examples
    ========

    >>> from sympy import Rational
    >>> from sympy.ntheory.egyptian_fraction import egypt_graham_jewett
    >>> egypt_graham_jewett(Rational(2,3))
    [3, 4, 12]
    >>> egypt_graham_jewett(Rational(3,7))
    [7, 8, 9, 56, 57, 72, 3192]

    See Also
    ========

    egypt_takenouchi(Rational) : Uses Takenouchi Algorithm
    egypt_greedy(Rational) : Uses Fibonacci-Sylvester Algorithm
    """

    x, y = int(Rational.p), int(Rational.q)
    l = [y] * x

    #l is now a list of integers whose reciprocals sum to x/y.
    #We shall now proceed to manipulate the elements of l without
    #changing the reciprocated sum until all elements are unique.

    while len(l) != len(set(l)):
        l.sort() #So the list has duplicates. Find a smallest pair
        for i in range(len(l) - 1):
            if l[i] == l[i + 1]:
                break
        #We have now identified a pair of identical elements: l[i] and l[i+1].
        #Now comes the application of the result of Graham and Jewett:
        l[i + 1] = l[i] + 1
        l.append(l[i]*(l[i] + 1)) #And we just iterate that until the list has no duplicates.  Ta da!
    return l

def egypt_takenouchi(Rational):
    """
    The algorithm suggested by Takenouchi (1921).
    Differs from the Graham-Jewett algorithm only in the handling of duplicates.

    References
    ==========

    - http://www.ics.uci.edu/~eppstein/numth/egypt/conflict.html

    Examples
    ========

    >>> from sympy import Rational
    >>> from sympy.ntheory.egyptian_fraction import egypt_takenouchi
    >>> egypt_takenouchi(Rational(3,7))
    [4, 28, 7]
    >>> egypt_takenouchi(Rational(7,23))
    [6, 12, 23, 138, 276]

    See Also
    ========

    egypt_greedy(Rational) : Uses Fibonacci-Sylvester Algorithm
    egypt_graham_jewett(Rational) : Uses Graham and Jewett Algorithm
    """

    x, y = int(Rational.p), int(Rational.q)
    l = [y] * x
    while len(l) != len(set(l)):
        l.sort()
        for i in range(len(l) - 1):
            if l[i] == l[i + 1]:
                break
        k = l[i]
        if k % 2 == 0:
            l[i] = (l[i] // 2)
            del l[i + 1]
        else:
            l[i],l[i + 1] = ((k + 1)//2), (k*(k + 1)//2)
    return l
