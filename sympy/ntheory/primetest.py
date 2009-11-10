"""
Primality testing

"""

_tiny_primes = [2, 3, 5, 7, 11, 13, 17, 19]
_max_tiny_prime = 19
_tiny_primes_set = set(_tiny_primes)

def _is_tiny_prime(n):
    return n <= _max_tiny_prime and n in _tiny_primes_set

def _has_tiny_factor(n):
    for p in _tiny_primes:
        if n % p == 0:
            return True
    return False

def _factor_pow2(n):
    """Factor powers of two from n. Return (s, t), with t odd, such
    that n = 2**s * t."""
    s, t = 0, n
    while not t & 1:
        t >>= 1
        s += 1
    return s, t

def _test(n, base):
    """Miller-Rabin strong pseudoprime test for one base.
    Return False if n is definitely composite, True if n is
    probably prime, with a probability greater than 3/4."""
    n = int(n)
    if n < 2:
        return False
    s, t = _factor_pow2(n-1)
    b = pow(base, t, n)
    if b == 1 or b == n-1:
        return True
    else:
        for j in xrange(1, s):
            b = (b**2) % n
            if b == n-1:
                return True
    return False

def mr(n, bases):
    """Perform a Miller-Rabin strong pseudoprime test on n using a
    given list of bases/witnesses.

    Reference:
    Richard Crandall & Carl Pomerance (2005), "Prime Numbers:
    A Computational Perspective", Springer, 2nd edition, 135-138
    """
    n = int(n)
    for base in bases:
        if not _test(n, base):
            return False
    return True

def mr_safe(n):
    """For n < 1e16, use the Miller-Rabin test to determine with
    certainty (unless the code is buggy!) whether n is prime.

    Reference for the bounds:
    http://primes.utm.edu/prove/prove2_3.html
    """
    n = int(n)
    if n < 1373653: return mr(n, [2, 3])
    if n < 25326001: return mr(n, [2, 3, 5])
    if n < 2152302898747: return mr(n, [2, 3, 5, 7, 11])
    if n < 3474749660383: return mr(n, [2, 3, 5, 7, 11, 13])
    if n < 341550071728321: return mr(n, [2, 3, 5, 7, 11, 13, 17])
    if n < 10000000000000000: return mr(n, [2, 3, 7, 61, 24251])
    raise ValueError("n too large")

def isprime(n):
    """
    Test whether n is a prime number. Negative primes (e.g. -2) are not
    considered prime. The function first looks for trivial factors,
    and if none is found, performs a Miller-Rabin strong pseudoprime
    test.

    Example usage
    =============
    >>> from sympy.ntheory import isprime
    >>> isprime(13)
        True
        >>> isprime(15)
        False

    """
    n = int(n)
    if n < 2: return False
    if n & 1 == 0: return n == 2
    if _is_tiny_prime(n):   return True
    if _has_tiny_factor(n): return False
    try:
        return mr_safe(n)
    except ValueError:
        # Should be good enough for practical purposes
        return mr(n, [2,3,5,7,11,13,17,19])
