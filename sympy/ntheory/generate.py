"""
Generating and counting primes.
"""

import random
from primetest import isprime


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
        # Use a binary search
        a = 0
        b = len(self._list)-1
        while 1:
            m = (a + b) // 2
            if self._list[m] == n:
                return m+1, m+1
            elif self._list[m] > n:
                b = m - 1
            elif self._list[m] < n:
                a = m + 1
            if a > b:
                return b+1, a+1

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
    """Return the nth prime, with the primes indexed as prime(1) = 2,
    prime(2) = 3, etc."""
    assert n > 0
    return sieve[n]

def primepi(n):
    """
    Return the value of the prime counting function pi(n) = the number
    of prime numbers less than or equal to n. The number n need not
    necessarily be an integer.
    """
    if n < 2:
        return 0
    else:
        n = int(n)
        return sieve.search(n)[0]


# TODO: the following functions could use the sieve for small n

def nextprime(n):
    """Return the smallest prime greater than n."""
    n = max(n, 0)
    while 1:
        n += 1
        if isprime(n):
            return n

def prevprime(n):
    """Return the largest prime smaller than n."""
    if n < 3:
        raise ValueError("no preceding primes")
    while 1:
        n -= 1
        if isprime(n):
            return n

def primerange(a, b):
    """Generate a list of all prime numbers in the range [a, b)."""
    assert a <= b
    a -= 1
    while 1:
        a = nextprime(a)
        if a < b:
            yield a
        else:
            return

def randprime(a, b):
    """Return a random prime number in the range [a, b)"""
    n = random.randint(a-1, b)
    p = nextprime(n)
    if p >= b:
        p = prevprime(b)
    if p < a:
        raise ValueError("no primes exist in the specified range")
    return p
