"""
Generating and counting primes.
"""

import random
from bisect import bisect
from primetest import isprime
from sympy.core.numbers import integer_nthroot


# Using arrays for sieving instead of lists greatly reduces
# memory consumption
from array import array as _array

def _arange(a, b):
    ar = _array('l', [0]*(b-a))
    for i, e in enumerate(xrange(a, b)):
        ar[i] = e
    return ar


class Sieve:
    """An infinite list of prime numbers, implemented as a dynamically
    growing sieve of Eratosthenes. When a lookup is requested involving
    a number that has not been sieved, the sieve is automatically
    extended up to that number."""

    _list = _array('l', [2, 3, 5, 7, 11, 13])

    def __repr__(self):
        return "<Sieve with %i primes sieved: 2, 3, 5, ... %i, %i>" % \
            (len(self._list), self._list[-2], self._list[-1])

    def extend(self, N):
        """Grow the sieve to cover all numbers <= N."""
        if N <= self._list[-1]:
            return

        # We need to sieve against all bases up to sqrt(n). If there
        # are too few, extend the list recursively.
        maxbase = int(N**0.5)+1
        self.extend(maxbase)

        # Create a new sieve starting from N**0.5
        begin = self._list[-1] + 1
        newsieve = _arange(begin, N+1)

        # Now eliminate all multiples of primes in [2, N**0.5]
        for p in self.primerange(2, maxbase):
            # Start counting at a multiple of p, offsetting
            # the index to account for the new sieve's base index
            startindex = (-begin) % p
            for i in xrange(startindex, len(newsieve), p):
                newsieve[i] = 0

        # Merge the sieves
        self._list += _array('l', [x for x in newsieve if x])

    def extend_to_no(self, n):
        """Extend to include (at least) the nth prime number"""
        while len(self._list) < n:
            self.extend(int(self._list[-1] * 1.5))

    def primerange(self, a, b):
        """Generate all prime numbers in the range [a, b)."""
        assert a <= b
        if b < 2:
            return
        a = max(2, a)
        self.extend(b)
        i = self.search(a)[1]
        maxi = len(self._list) + 1
        while i < maxi:
            p = self._list[i-1]
            if p < b:
                yield p
                i += 1
            else:
                return

    def search(self, n):
        """For n >= 2, return the tightest a, b such that
        self[a] <= n <= self[b]"""
        assert n >= 2
        if n > self._list[-1]:
            self.extend(n)
        b = bisect(self._list, n)
        if self._list[b-1] == n:
            return b, b
        else:
            return b, b+1

    def __contains__(self, n):
        if n < 2:
            return False
        a, b = self.search(n)
        return a == b

    def __getitem__(self, n):
        """Return the nth prime number"""
        self.extend_to_no(n)
        return self._list[n-1]


# Generate a global object for repeated use in trial division etc
sieve = Sieve()


def prime(n):
    """ Return the nth prime, with the primes indexed as prime(1) = 2,
        prime(2) = 3, etc.... The nth prime is approximately n*log(n) and
        can never be larger than 2**n.

        Reference: http://primes.utm.edu/glossary/xpage/BertrandsPostulate.html
    """

    assert n > 0
    return sieve[n]

def primepi(n):
    """ Return the value of the prime counting function pi(n) = the number
        of prime numbers less than or equal to n. The number n need not
        necessarily be an integer.
    """

    if n < 2:
        return 0
    else:
        n = int(n)
        return sieve.search(n)[0]

def nextprime(n, i=1):
    """ Return the ith prime greater than n.

        Potential primes are located at 6*j +/- 1.

        >>> from sympy import nextprime
        >>> [(i, nextprime(i)) for i in range(10, 15)]
        [(10, 11), (11, 13), (12, 13), (13, 17), (14, 17)]
        >>> nextprime(2, i=2) # the 2nd prime after 2
        5
    """

    if i > 1:
        pr = n
        j = 1
        while 1:
            pr = nextprime(pr)
            j += 1
            if j > i:
                break
        return pr

    n = int(n)
    if n < 2:
        return 2
    if n < 7:
        return {2: 3, 3: 5, 4: 5, 5: 7, 6: 7}[n]
    nn = 6*(n//6)
    if nn == n:
        n += 1
        if isprime(n):
            return n
        n += 4
    elif n - nn == 5:
        n += 2
        if isprime(n):
            return n
        n += 4
    else:
        n = nn + 5
    while 1:
        if isprime(n):
            return n
        n += 2
        if isprime(n):
            return n
        n += 4

def prevprime(n):
    """ Return the largest prime smaller than n.

        Potential primes are located at 6*j +/- 1.

        >>> from sympy import prevprime
        >>> [(i, prevprime(i)) for i in range(10, 15)]
        [(10, 7), (11, 7), (12, 11), (13, 11), (14, 13)]
    """

    n = int(n)
    if n < 3:
        raise ValueError("no preceding primes")
    if n < 8:
        return {3: 2, 4: 3, 5: 3, 6: 5, 7: 5}[n]
    nn = 6*(n//6)
    if n - nn <= 1:
        n = nn - 1
        if isprime(n):
            return n
        n -= 4
    else:
        n = nn + 1
    while 1:
        if isprime(n):
            return n
        n -= 2
        if isprime(n):
            return n
        n -= 4

def primerange(a, b):
    """ Generate a list of all prime numbers in the range [a, b).

        Some famous conjectures about the occurence of primes in a given
        range are [1]:

        - Twin primes: though often not, the following will give 2 primes an oo
                      number of times:
                      primerange(6*n - 1, 6*n + 2)
        - Legendre's: the following always yields at least one prime
                      primerange(n**2, (n+1)**2+1)
        - Bertrand's (proven): there is always a prime in the range
                      primerange(n, 2*n)
        - Brocard's: there are at least four primes in the range
                     primerange(prime(n)**2, prime(n+1)**2)

        The average gap between primes is log(n) [2];
        the gap between primes can be arbitrarily large since sequences of
        composite numbers are arbitrarily large, e.g. the numbers in the sequence
        n!+2, n!+3 ... n!+n are all composite.

        References:
            [1] http://en.wikipedia.org/wiki/Prime_number
            [2] http://primes.utm.edu/notes/gaps.html
    """
    assert a <= b
    a -= 1
    while 1:
        a = nextprime(a)
        if a < b:
            yield a
        else:
            return

def randprime(a, b):
    """ Return a random prime number in the range [a, b).

        Bertrand's postulate assures that
        randprime(a, 2*a) will always succeed for a > 1.

        Reference: http://en.wikipedia.org/wiki/Bertrand's_postulate
    """

    n = random.randint(a-1, b)
    p = nextprime(n)
    if p >= b:
        p = prevprime(b)
    if p < a:
        raise ValueError("no primes exist in the specified range")
    return p

def primorial(n, nth=True):
    """ Returns the product of either a) the first n primes (default) or
        b) the primes less than or equal to n (when `nth`=False).

        >>> from sympy.ntheory.generate import primorial, randprime, primerange
        >>> from sympy import factorint, Mul, primefactors
        >>> primorial(4) # the first 4 primes are 2, 3, 5, 7
        210
        >>> primorial(4, nth=0) # primes <= 4 are 2 and 3
        6
        >>> primorial(1)
        2
        >>> primorial(1, nth=0)
        1

        One can argue that the primes are infinite since if you take
        a set of primes and multiply them together (e.g. the primorial) and
        then add or subtract 1, the result cannot be divided by any of the
        original factors, hence either 1 or more primes must divide this
        product of primes.

        >>> factorint(primorial(4) + 1)
        {211: 1}
        >>> factorint(primorial(4) - 1)
        {11: 1, 19: 1}
        >>> p = list(primerange(10, 20))
        >>> sorted(set(primefactors(Mul(*p) + 1)).difference(set(p)))
        [2, 5, 31, 149]
    """

    if n < 1:
        raise ValueError("primorial argument must be >= 1")
    p = 1
    if nth:
        for i in range(1, n + 1):
            p *= prime(i)
    else:
        for i in primerange(2, n + 1):
            p *= i
    return p

def cycle_length(f, x0, nmax=None, values=False):
    """For a given iterated sequence, return a generator that gives
    the length of the iterated cycle (lambda) and the length of terms
    before the cycle begins (mu); if ``values`` is True then the
    terms of the sequence will be returned instead.

    Note: more than the first lambda + mu terms may be returned and this
    is the cost of cycle detection with Brent's method; there are, however,
    generally less terms calculated than would have been calculated if the
    proper ending point were determined, e.g. by using Floyd's method.

    >>> from sympy.ntheory.generate import cycle_length
    >>> from random import Random

    This will yield successive values of i <-- func(i):

        >>> def iter(func, i):
        ...     while 1:
        ...         ii = func(i)
        ...         yield ii
        ...         i = ii
        ...

    A function is defined:

        >>> func = lambda i: (i**2 + 1) % 51

    and given a seed of 2 and the mu and lambda terms calculated:

        >>> cycle_length(func, 4).next()
        (6, 2)

    We can see what is meant by looking at the output:

        >>> n = cycle_length(func, 4, values=True)
        >>> list(ni for ni in n)
        [17, 35, 2, 5, 26, 14, 44, 50, 2, 5, 26, 14]

                  \_______________/
                   6 values after
                    the first 2

    If a sequence is suspected of being longer than you might wish, ``nmax``
    can be used to exit early (in which mu will be returned as None:

        >>> cycle_length(func, 4, nmax = 4).next()
        (4, None)
        >>> [ni for ni in cycle_length(func, 4, nmax = 4, values=True)]
        [17, 35, 2, 5]

    Code modified from:
        http://en.wikipedia.org/wiki/Cycle_detection.
    """

    # main phase: search successive powers of two
    power = lam = 1
    tortoise, hare = x0, f(x0) # f(x0) is the element/node next to x0.
    i = 0
    while tortoise != hare and (not nmax or i < nmax):
        i += 1
        if power == lam:   # time to start a new power of two?
            tortoise = hare
            power *= 2
            lam = 0
        if values:
            yield hare
        hare = f(hare)
        lam += 1
    if nmax and i == nmax:
        if values:
            return
        else:
            yield nmax, None
            return
    if not values:
        # Find the position of the first repetition of length lambda
        mu = 0
        tortoise = hare = x0
        for i in range(lam):
            hare = f(hare)
        while tortoise != hare:
            tortoise = f(tortoise)
            hare = f(hare)
            mu += 1
        if mu:
            mu -= 1
        yield lam, mu
