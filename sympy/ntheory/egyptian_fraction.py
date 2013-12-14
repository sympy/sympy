from __future__ import print_function, division

from sympy import Rational
from fractions import gcd
from random import randint

def egyptian(rat, choice):
    """
    Return the denominators of an Egtyptian fraction expansion of
    the said rational.

    The egypt function takes a Rational type and choice as input. 
    The choice defines which algorithm to be used to generate the
    denominators. Currently the following algorithms are 
    supported :
    
    1) Greedy Algorithm :
    Greedy algorithm for Egyptian fraction expansion
    Also called the Fibonacci-Sylvester algorithm
    At each step, extract the largest unit fraction less
    than the target and replace the target with the remainder

    2) Graham Jewett Algorithm :
    The algorithm suggested by the result of Graham and Jewett.
    Note that this has a tendency to blow up: the length of the 
    resulting expansion is always 2**(x/gcd(x, y)) - 1

    3) Takenouchi Algorithm :
    The algorithm suggested by Takenouchi (1921).
    Differs from the Graham-Jewett algorithm only in the handling 
    of duplicates.

    References
    =========
    
    - https://en.wikipedia.org/wiki/Greedy_algorithm_for_Egyptian_fractions
    - http://en.wikipedia.org/wiki/Egyptian_fraction#Modern_number_theory
    - http://www.ics.uci.edu/~eppstein/numth/egypt/conflict.html
    
    Examples
    ========

    >>> from sympy import Rational
    >>> from sympy.ntheory.egyptian_fraction import egyptian
    >>> egyptian(Ratioanl(3, 7), "Greedy")
    [3, 11, 231]
    >>> egyptian(Rational(3, 7), "Graham Jewett")
    [7, 8, 9, 56, 57, 72, 3192]
    >>> egyptian(Rational(3, 7), "Takenouchi")
    [4, 28, 7]

    """
    def egypt_greedy():

        x, y = rat.as_numer_denom()
        if x == 1:
            return [y]
        else:
            a = (-y) % (x)
            b = y*((y//x) + 1)
            c = gcd(a, b)
            if c > 1:
                num, denom = (a//c), (b//c)
            else:
                num, denom = a, b
            return [(y//x) + 1] + egyptian(Rational(num, denom), "Greedy")

    def egypt_graham_jewett():

        x, y = rat.as_numer_denom()
        l = [y] * x

        # l is now a list of integers whose reciprocals sum to x/y.
        # We shall now proceed to manipulate the elements of l without
        # changing the reciprocated sum until all elements are unique.

        while len(l) != len(set(l)):
            l.sort() # So the list has duplicates. Find a smallest pair
            for i in range(len(l) - 1):
                if l[i] == l[i + 1]:
                    break
            # We have now identified a pair of identical 
            # elements: l[i] and l[i + 1].
            # Now comes the application of the result of Graham and Jewett:
            l[i + 1] = l[i] + 1 
            # And we just iterate that until the list has no duplicates.
            l.append(l[i]*(l[i] + 1))
        return l

    def egypt_takenouchi():

        x, y = rat.as_numer_denom()
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
                l[i], l[i + 1] = ((k + 1)//2), (k*(k + 1)//2)
        return l
    
    if choice == "Greedy":
        return egypt_greedy()
    elif choice == "Graham Jewett":
        return egypt_graham_jewett()
    elif choice == "Takenouchi":
        return egypt_takenouchi()
    else :
        raise ValueError("Entered Invalid Algorithm")
