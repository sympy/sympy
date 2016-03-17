"""
Generating and counting primes.

"""
from __future__ import print_function, division

import random
from bisect import bisect
# Using arrays for sieving instead of lists greatly reduces
# memory consumption
from array import array as _array

from .primetest import isprime
from sympy.core.compatibility import as_int, range


def _arange(a, b):
    ar = _array('l', [0]*(b - a))
    for i, e in enumerate(range(a, b)):
        ar[i] = e
    return ar


class Sieve:
    """An infinite list of prime numbers, implemented as a dynamically
    growing sieve of Eratosthenes. When a lookup is requested involving
    an odd number that has not been sieved, the sieve is automatically
    extended up to that number.

    >>> from sympy import sieve
    >>> from array import array # this line and next for doctest only
    >>> sieve._list = array('l', [2, 3, 5, 7, 11, 13])

    >>> 25 in sieve
    False
    >>> sieve._list
    array('l', [2, 3, 5, 7, 11, 13, 17, 19, 23])
    """

    # data shared (and updated) by all Sieve instances
    _list = _array('l', [2, 3, 5, 7, 11, 13])

    def __repr__(self):
        return "<Sieve with %i primes sieved: 2, 3, 5, ... %i, %i>" % \
            (len(self._list), self._list[-2], self._list[-1])

    def extend(self, n):
        """Grow the sieve to cover all primes <= n (a real number).

        Examples
        ========

        >>> from sympy import sieve
        >>> from array import array # this line and next for doctest only
        >>> sieve._list = array('l', [2, 3, 5, 7, 11, 13])

        >>> sieve.extend(30)
        >>> sieve[10] == 29
        True
        """
        n = int(n)
        if n <= self._list[-1]:
            return

        # We need to sieve against all bases up to sqrt(n).
        # This is a recursive call that will do nothing if there are enough
        # known bases already.
        maxbase = int(n**0.5) + 1
        self.extend(maxbase)

        # Create a new sieve starting from sqrt(n)
        begin = self._list[-1] + 1
        newsieve = _arange(begin, n + 1)

        # Now eliminate all multiples of primes in [2, sqrt(n)]
        for p in self.primerange(2, maxbase):
            # Start counting at a multiple of p, offsetting
            # the index to account for the new sieve's base index
            startindex = (-begin) % p
            for i in range(startindex, len(newsieve), p):
                newsieve[i] = 0

        # Merge the sieves
        self._list += _array('l', [x for x in newsieve if x])

    def extend_to_no(self, i):
        """Extend to include the ith prime number.

        i must be an integer.

        The list is extended by 50% if it is too short, so it is
        likely that it will be longer than requested.

        Examples
        ========

        >>> from sympy import sieve
        >>> from array import array # this line and next for doctest only
        >>> sieve._list = array('l', [2, 3, 5, 7, 11, 13])

        >>> sieve.extend_to_no(9)
        >>> sieve._list
        array('l', [2, 3, 5, 7, 11, 13, 17, 19, 23])
        """
        i = as_int(i)
        while len(self._list) < i:
            self.extend(int(self._list[-1] * 1.5))

    def primerange(self, a, b):
        """Generate all prime numbers in the range [a, b).

        Examples
        ========

        >>> from sympy import sieve
        >>> print([i for i in sieve.primerange(7, 18)])
        [7, 11, 13, 17]
        """
        from sympy.functions.elementary.integers import ceiling

        # wrapping ceiling in int will raise an error if there was a problem
        # determining whether the expression was exactly an integer or not
        a = max(2, int(ceiling(a)))
        b = int(ceiling(b))
        if a >= b:
            return
        self.extend(b)
        i = self.search(a)[1]
        maxi = len(self._list) + 1
        while i < maxi:
            p = self._list[i - 1]
            if p < b:
                yield p
                i += 1
            else:
                return

    def search(self, n):
        """Return the indices i, j of the primes that bound n.

        If n is prime then i == j.

        Although n can be an expression, if ceiling cannot convert
        it to an integer then an n error will be raised.

        Examples
        ========

        >>> from sympy import sieve
        >>> sieve.search(25)
        (9, 10)
        >>> sieve.search(23)
        (9, 9)
        """
        from sympy.functions.elementary.integers import ceiling

        # wrapping ceiling in int will raise an error if there was a problem
        # determining whether the expression was exactly an integer or not
        test = int(ceiling(n))
        n = int(n)
        if n < 2:
            raise ValueError("n should be >= 2 but got: %s" % n)
        if n > self._list[-1]:
            self.extend(n)
        b = bisect(self._list, n)
        if self._list[b - 1] == test:
            return b, b
        else:
            return b, b + 1

    def __contains__(self, n):
        try:
            n = as_int(n)
            assert n >= 2
        except (ValueError, AssertionError):
            return False
        if n % 2 == 0:
            return n == 2
        a, b = self.search(n)
        return a == b

    def __getitem__(self, n):
        """Return the nth prime number"""
        if isinstance(n, slice):
            self.extend_to_no(n.stop)
            return self._list[n.start - 1:n.stop - 1:n.step]
        else:
            n = as_int(n)
            self.extend_to_no(n)
            return self._list[n - 1]

# Generate a global object for repeated use in trial division etc
sieve = Sieve()


def prime(nth):
    """ Return the nth prime, with the primes indexed as prime(1) = 2,
        prime(2) = 3, etc.... The nth prime is approximately n*log(n) and
        can never be larger than 2**n.

        References
        ==========

        - http://primes.utm.edu/glossary/xpage/BertrandsPostulate.html

        Examples
        ========

        >>> from sympy import prime
        >>> prime(10)
        29
        >>> prime(1)
        2

        See Also
        ========

        sympy.ntheory.primetest.isprime : Test if n is prime
        primerange : Generate all primes in a given range
        primepi : Return the number of primes less than or equal to n
    """
    n = as_int(nth)
    if n < 1:
        raise ValueError("nth must be a positive integer; prime(1) == 2")
    return sieve[n]


def primepi(n):
    """ Return the value of the prime counting function pi(n) = the number
        of prime numbers less than or equal to n.

        Examples
        ========

        >>> from sympy import primepi
        >>> primepi(25)
        9
        >>> primepi(100)
        25
        >>> primepi(10000000000)
        455052511
        See Also
        ========

        sympy.ntheory.primetest.isprime : Test if n is prime
        primerange : Generate all primes in a given range
        prime : Return the nth prime
    """
    n = int(n)
    if n<2:
    	return 0	
    lim = int(n**0.5) #lim = floor(sqrt(n))
    lim-=1
    lim = max(lim,0)
    while lim*lim<=n:
        lim+=1
    lim-=1
    arr1 = [0]*(lim+1) # stores values of pi(i)
    arr2 = [0]*(lim+1) # stores values of pi(n/i)
    for i in range(1,lim+1):
    	arr1[i] = i - 1
    	if i == 0:
    	    arr2[i] = 0
    	else:
    	    arr2[i] = n // i -1	
    for i in range(2,lim+1):
    	if arr1[i] == arr1[i-1]:
    	    continue
    	p = arr1[i-1]
    	for j in range(1,min(n//(i*i),lim)+1):
    	    st = i*j
    	    if st<=lim:
    	        arr2[j]-=arr2[st]-p	
    	    else:
    	    	arr2[j]-=arr1[n//st]-p
        for j in range(lim,i*i-1,-1):
            arr1[j]-=arr1[j//i]-p	
    return arr2[1] #pi(n)


def nextprime(n, ith=1):
    """ Return the ith prime greater than n.

        i must be an integer.

        Notes
        =====

        Potential primes are located at 6*j +/- 1. This
        property is used during searching.

        >>> from sympy import nextprime
        >>> [(i, nextprime(i)) for i in range(10, 15)]
        [(10, 11), (11, 13), (12, 13), (13, 17), (14, 17)]
        >>> nextprime(2, ith=2) # the 2nd prime after 2
        5

        See Also
        ========

        prevprime : Return the largest prime smaller than n
        primerange : Generate all primes in a given range

    """
    n = int(n)
    i = as_int(ith)
    if i > 1:
        pr = n
        j = 1
        while 1:
            pr = nextprime(pr)
            j += 1
            if j > i:
                break
        return pr

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

        Notes
        =====

        Potential primes are located at 6*j +/- 1. This
        property is used during searching.

        >>> from sympy import prevprime
        >>> [(i, prevprime(i)) for i in range(10, 15)]
        [(10, 7), (11, 7), (12, 11), (13, 11), (14, 13)]

        See Also
        ========

        nextprime : Return the ith prime greater than n
        primerange : Generates all primes in a given range
    """
    from sympy.functions.elementary.integers import ceiling

    # wrapping ceiling in int will raise an error if there was a problem
    # determining whether the expression was exactly an integer or not
    n = int(ceiling(n))
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

        If the range exists in the default sieve, the values will
        be returned from there; otherwise values will be returned
        but will not modify the sieve.

        Notes
        =====

        Some famous conjectures about the occurence of primes in a given
        range are [1]:

        - Twin primes: though often not, the following will give 2 primes
                    an infinite number of times:
                        primerange(6*n - 1, 6*n + 2)
        - Legendre's: the following always yields at least one prime
                        primerange(n**2, (n+1)**2+1)
        - Bertrand's (proven): there is always a prime in the range
                        primerange(n, 2*n)
        - Brocard's: there are at least four primes in the range
                        primerange(prime(n)**2, prime(n+1)**2)

        The average gap between primes is log(n) [2]; the gap between
        primes can be arbitrarily large since sequences of composite
        numbers are arbitrarily large, e.g. the numbers in the sequence
        n! + 2, n! + 3 ... n! + n are all composite.

        References
        ==========

        1. http://en.wikipedia.org/wiki/Prime_number
        2. http://primes.utm.edu/notes/gaps.html

        Examples
        ========

        >>> from sympy import primerange, sieve
        >>> print([i for i in primerange(1, 30)])
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

        The Sieve method, primerange, is generally faster but it will
        occupy more memory as the sieve stores values. The default
        instance of Sieve, named sieve, can be used:

        >>> list(sieve.primerange(1, 30))
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

        See Also
        ========

        nextprime : Return the ith prime greater than n
        prevprime : Return the largest prime smaller than n
        randprime : Returns a random prime in a given range
        primorial : Returns the product of primes based on condition
        Sieve.primerange : return range from already computed primes
                           or extend the sieve to contain the requested
                           range.
    """
    from sympy.functions.elementary.integers import ceiling

    # if we already have the range, return it
    if b <= sieve._list[-1]:
        for i in sieve.primerange(a, b):
            yield i
        return
    # otherwise compute, without storing, the desired range
    if a >= b:
        return
    # wrapping ceiling in int will raise an error if there was a problem
    # determining whether the expression was exactly an integer or not
    a = int(ceiling(a)) - 1
    b = int(ceiling(b))
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

        References
        ==========

        - http://en.wikipedia.org/wiki/Bertrand's_postulate

        Examples
        ========

        >>> from sympy import randprime, isprime
        >>> randprime(1, 30) #doctest: +SKIP
        13
        >>> isprime(randprime(1, 30))
        True

        See Also
        ========

        primerange : Generate all primes in a given range

    """
    if a >= b:
        return
    a, b = map(int, (a, b))
    n = random.randint(a - 1, b)
    p = nextprime(n)
    if p >= b:
        p = prevprime(b)
    if p < a:
        raise ValueError("no primes exist in the specified range")
    return p


def primorial(n, nth=True):
    """
    Returns the product of the first n primes (default) or
    the primes less than or equal to n (when ``nth=False``).

    >>> from sympy.ntheory.generate import primorial, randprime, primerange
    >>> from sympy import factorint, Mul, primefactors, sqrt
    >>> primorial(4) # the first 4 primes are 2, 3, 5, 7
    210
    >>> primorial(4, nth=False) # primes <= 4 are 2 and 3
    6
    >>> primorial(1)
    2
    >>> primorial(1, nth=False)
    1
    >>> primorial(sqrt(101), nth=False)
    210

    One can argue that the primes are infinite since if you take
    a set of primes and multiply them together (e.g. the primorial) and
    then add or subtract 1, the result cannot be divided by any of the
    original factors, hence either 1 or more new primes must divide this
    product of primes.

    In this case, the number itself is a new prime:

    >>> factorint(primorial(4) + 1)
    {211: 1}

    In this case two new primes are the factors:

    >>> factorint(primorial(4) - 1)
    {11: 1, 19: 1}

    Here, some primes smaller and larger than the primes multiplied together
    are obtained:

    >>> p = list(primerange(10, 20))
    >>> sorted(set(primefactors(Mul(*p) + 1)).difference(set(p)))
    [2, 5, 31, 149]

    See Also
    ========

    primerange : Generate all primes in a given range

    """
    if nth:
        n = as_int(n)
    else:
        n = int(n)
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
    terms of the sequence will be returned instead. The sequence is
    started with value ``x0``.

    Note: more than the first lambda + mu terms may be returned and this
    is the cost of cycle detection with Brent's method; there are, however,
    generally less terms calculated than would have been calculated if the
    proper ending point were determined, e.g. by using Floyd's method.

    >>> from sympy.ntheory.generate import cycle_length

    This will yield successive values of i <-- func(i):

        >>> def iter(func, i):
        ...     while 1:
        ...         ii = func(i)
        ...         yield ii
        ...         i = ii
        ...

    A function is defined:

        >>> func = lambda i: (i**2 + 1) % 51

    and given a seed of 4 and the mu and lambda terms calculated:

        >>> next(cycle_length(func, 4))
        (6, 2)

    We can see what is meant by looking at the output:

        >>> n = cycle_length(func, 4, values=True)
        >>> list(ni for ni in n)
        [17, 35, 2, 5, 26, 14, 44, 50, 2, 5, 26, 14]

    There are 6 repeating values after the first 2.

    If a sequence is suspected of being longer than you might wish, ``nmax``
    can be used to exit early (and mu will be returned as None):

        >>> next(cycle_length(func, 4, nmax = 4))
        (4, None)
        >>> [ni for ni in cycle_length(func, 4, nmax = 4, values=True)]
        [17, 35, 2, 5]

    Code modified from:
        http://en.wikipedia.org/wiki/Cycle_detection.
    """

    nmax = int(nmax or 0)

    # main phase: search successive powers of two
    power = lam = 1
    tortoise, hare = x0, f(x0)  # f(x0) is the element/node next to x0.
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
