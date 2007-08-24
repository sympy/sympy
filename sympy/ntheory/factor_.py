"""
Integer factorization
"""

from sympy.core import numbers
import random
from primetest import isprime
from generate import sieve, prime


def multiplicity(p, n):
    """
    Return the multiplicity of the number p in n; that is, the greatest
    number m such that p**m divides n.

    Example usage
    =============
        >>> multiplicity(5, 8)
        0
        >>> multiplicity(5, 5)
        1
        >>> multiplicity(5, 25)
        2
        >>> multiplicity(5, 125)
        3
        >>> multiplicity(5, 250)
        3
    """
    m = 0
    quot, rem = divmod(n, p)
    while rem == 0:
        quot, rem = divmod(quot, p)
        m += 1
    return m


# Currently unused by factorint()
def pollard_rho(n, max_iters=5, seed=1234):
    """Use Pollard's rho method to try to extract a factor of n. The
    returned factor may be a composite number. A maximum of max_iters
    iterations are performed; if no factor is found, None is returned.

    The rho algorithm is a Monte Carlo method whose outcome can
    be affected by changing the random seed value.

    References
    ==========
    Richard Crandall & Carl Pomerance (2005), "Prime Numbers:
    A Computational Perspective", Springer, 2nd edition, 229-231

    """
    random.seed(seed + max_iters)
    for i in range(max_iters):
        a = random.randint(1, n-3)
        s = random.randint(0, n-1)
        U = V = s
        F = lambda x: (x**2 + a) % a
        while 1:
            U = F(U)
            V = F(F(V))
            g = numbers.gcd(U-V, n)
            if g == 1:
                continue
            if g == n:
                break
            return g
    return None


def pollard_pm1(n, B=10, seed=1234):
    """Use Pollard's p-1 method to try to extract a factor of n. The
    returned factor may be a composite number. The search is performed
    up to a smoothness bound B; if no factor is found, None is
    returned.

    The p-1 algorithm is a Monte Carlo method whose outcome can
    be affected by changing the random seed value.

    Example usage
    =============
    With the default smoothness bound, this number can't be cracked:
        >>> pollard_pm1(21477639576571)

    Increasing the smoothness bound helps:
        >>> pollard_pm1(21477639576571, 2000)
        4410317L

    References
    ==========
    Richard Crandall & Carl Pomerance (2005), "Prime Numbers:
    A Computational Perspective", Springer, 2nd edition, 236-238

    """
    from math import log
    random.seed(seed + B)
    a = random.randint(2, n-1)
    for p in sieve.primerange(2, B):
        e = int(log(B, p))
        a = pow(a, p**e, n)
    g = numbers.gcd(a-1, n)
    if 1 < g < n:
        return g
    else:
        return None


def trial(n, candidates=None):
    """
    Factor n as far as possible through trial division, taking
    candidate factors from the given list. If no list of candidate
    factors is given, the prime numbers in the interval [2, sqrt(n)]
    are used, which guarantees a complete factorization.

    The returned value is a list [(p1, e1), ...] such that
    n = p1**e1 * p2**e2 * ... If n could not be completely factored
    using numbers in the given range, the last p might be composite.

    Example usage
    =============

    A complete factorization:

        >>> trial(36960)
        [(2, 5), (3, 1), (5, 1), (7, 1), (11, 1)]

    This won't find the factors 7 and 11:

        >>> trial(36960, [2, 3, 5])
        [(2, 5), (3, 1), (5, 1), (77, 1)]

    """
    if n == 1:
        return []
    if candidates is None:
        candidates = sieve.primerange(2, int(n**0.5)+1)
    factors = []
    for k in candidates:
        m = multiplicity(k, n)
        if m != 0:
            n //= k**m
            factors = factors + [(k, m)]
        if isprime(n):
            return factors + [(int(n), 1)]
        elif n == 1:
            return factors
    return factors + [(int(n), 1)]


def factorint(n, limit=None, verbose=False):
    """
    Given a positive integer n, factorint(n) returns a list
    [(p_1, m_1), (p_2, m_2), ...] with all p prime and n = p_1**m_1 *
    p_2**m_2 * ...

    Special cases: 1 factors as [], 0 factors as [(0, 1)], and negative
    integers factor as [(-1, 1), ...].

    The function uses a composite algorithm, switching between
    Pollard's p-1 method and looking for small factors through trial
    division.

    It is sometimes useful to look only for small factors. If 'limit'
    is specified, factorint will only perform trial division with
    candidate factors up to this limit (and p-1 search up to the same
    smoothness bound). As a result, the last 'prime' in the returned
    list may be composite.

    Example usage
    =============

    Here are some simple factorizations (with at most six digits in the
    second largest factor). They should all complete within a fraction
    of a second:

        >>> factorint(1)
        []

        >>> factorint(100)
        [(2, 2), (5, 2)]

        >>> factorint(17*19)
        [(17, 1), (19, 1)]

        >>> factorint(prime(100)*prime(1000)*prime(10000))
        [(541, 1), (7919, 1), (104729, 1)]

        >>> factorint(2**(2**6) + 1)
        [(274177, 1), (67280421310721L, 1)]

    Factors on the order of 10 digits can generally be found quickly.
    The following computations should complete within a few seconds:

        >>> factorint(21477639576571)
        [(4410317, 1), (4869863, 1)]

        >>> factorint(12345678910111213141516)
        [(2, 2), (2507191691L, 1), (1231026625769L, 1)]

        >>> factorint(5715365922033905625269)
        [(74358036521L, 1), (76862786989L, 1)]

    This number has an enormous semiprime factor that is better
    ignored:

        >>> a = 1407633717262338957430697921446883
        >>> factorint(a, limit=10000)
        [(7, 1), (991, 1), (202916782076162456022877024859L, 1)]
        >>> isprime(_[-1][0])
        False

    """
    n = int(n)

    if n < 0: return [(-1, 1)] + factorint(-n, limit)
    if n == 0:
        return [(0, 1)]
    if n == 1:
        return []
    if isprime(n):
        return [(n, 1)]
    if limit is None:
        limit = int(n**0.5) + 1

    factors = []
    low, high = 2, 50

    while 1:
        # Trial divide for small factors first
        tfactors = trial(n, sieve.primerange(low, min(high, limit)))

        if verbose:
            print "trial division from", low, "to", \
                min(high,limit)-1, "gave", tfactors

        # If all were primes, we're done
        if isprime(tfactors[-1][0]):
            factors += tfactors
            break

        elif tfactors[-1][0] == 1:
            factors += tfactors[:-1]
            break
        else:
            factors += tfactors[:-1]
            n = tfactors[-1][0]

        # If we're lucky, Pollard's p-1 will extract a large factor
        w = pollard_pm1(n, high)
        if verbose:
            print "pollard p-1 with smoothness bound", high, "gave", w
            print

        if w is not None:
            # w may be composite
            for f, m in factorint(w, limit):
                m *= multiplicity(f, n)
                factors += [(f, m)]
                n //= f**(m)

        if n == 1:
            break

        if isprime(n):
            factors += [(int(n), 1)]
            break

        if high > limit:
            factors += [(int(n), 1)]
            break

        low, high = high, high*5

    return sorted(factors)


def primefactors(n, limit=None, verbose=False):
    """Return a list of n's prime factors, ignoring multiplicity.
    Unlike factorint(), primefactors() only returns prime numbers;
    i.e., it does not return -1 or 0, and if 'limit' is set too
    low for all factors to be found, composite factors are ignored.

    Example usage
    =============

        >>> primefactors(6)
        [2, 3]
        >>> primefactors(-5)
        [5]

        >>> factorint(123456)
        [(2, 6), (3, 1), (643, 1)]
        >>> primefactors(123456)
        [2, 3, 643]

        >>> factorint(10000000001, limit=1000)
        [(101, 1), (99009901, 1)]
        >>> isprime(99009901)
        False
        >>> primefactors(10000000001, limit=1000)
        [101]

    """
    n = int(n)
    s = []
    factors = factorint(n, limit, verbose)
    for p, _ in factors[:-1]:
        if p not in [-1, 0, 1]:
            s += [p]
    if isprime(factors[-1][0]):
        s += [factors[-1][0]]
    return s
