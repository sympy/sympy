"""
Integer factorization
"""

from collections import defaultdict
from functools import reduce
import random
import math

from sympy.core import sympify
from sympy.core.containers import Dict
from sympy.core.evalf import bitcount
from sympy.core.expr import Expr
from sympy.core.function import Function
from sympy.core.logic import fuzzy_and
from sympy.core.mul import Mul
from sympy.core.numbers import igcd, ilcm, Rational, Integer
from sympy.core.power import integer_nthroot, Pow, integer_log
from sympy.core.singleton import S
from sympy.external.gmpy import SYMPY_INTS
from .primetest import isprime
from .generate import sieve, primerange, nextprime
from .digits import digits
from sympy.utilities.iterables import flatten
from sympy.utilities.misc import as_int, filldedent
from .ecm import _ecm_one_factor

# Note: This list should be updated whenever new Mersenne primes are found.
# Refer: https://www.mersenne.org/
MERSENNE_PRIME_EXPONENTS = (2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203,
 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049,
 216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011, 24036583,
 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933)

# compute more when needed for i in Mersenne prime exponents
PERFECT = [6]  # 2**(i-1)*(2**i-1)
MERSENNES = [3]  # 2**i - 1


def _ismersenneprime(n):
    global MERSENNES
    j = len(MERSENNES)
    while n > MERSENNES[-1] and j < len(MERSENNE_PRIME_EXPONENTS):
        # conservatively grow the list
        MERSENNES.append(2**MERSENNE_PRIME_EXPONENTS[j] - 1)
        j += 1
    return n in MERSENNES


def _isperfect(n):
    global PERFECT
    if n % 2 == 0:
        j = len(PERFECT)
        while n > PERFECT[-1] and j < len(MERSENNE_PRIME_EXPONENTS):
            # conservatively grow the list
            t = 2**(MERSENNE_PRIME_EXPONENTS[j] - 1)
            PERFECT.append(t*(2*t - 1))
            j += 1
    return n in PERFECT


small_trailing = [0] * 256
for j in range(1,8):
    small_trailing[1<<j::1<<(j+1)] = [j] * (1<<(7-j))


def smoothness(n):
    """
    Return the B-smooth and B-power smooth values of n.

    The smoothness of n is the largest prime factor of n; the power-
    smoothness is the largest divisor raised to its multiplicity.

    Examples
    ========

    >>> from sympy.ntheory.factor_ import smoothness
    >>> smoothness(2**7*3**2)
    (3, 128)
    >>> smoothness(2**4*13)
    (13, 16)
    >>> smoothness(2)
    (2, 2)

    See Also
    ========

    factorint, smoothness_p
    """

    if n == 1:
        return (1, 1)  # not prime, but otherwise this causes headaches
    facs = factorint(n)
    return max(facs), max(m**facs[m] for m in facs)


def smoothness_p(n, m=-1, power=0, visual=None):
    """
    Return a list of [m, (p, (M, sm(p + m), psm(p + m)))...]
    where:

    1. p**M is the base-p divisor of n
    2. sm(p + m) is the smoothness of p + m (m = -1 by default)
    3. psm(p + m) is the power smoothness of p + m

    The list is sorted according to smoothness (default) or by power smoothness
    if power=1.

    The smoothness of the numbers to the left (m = -1) or right (m = 1) of a
    factor govern the results that are obtained from the p +/- 1 type factoring
    methods.

        >>> from sympy.ntheory.factor_ import smoothness_p, factorint
        >>> smoothness_p(10431, m=1)
        (1, [(3, (2, 2, 4)), (19, (1, 5, 5)), (61, (1, 31, 31))])
        >>> smoothness_p(10431)
        (-1, [(3, (2, 2, 2)), (19, (1, 3, 9)), (61, (1, 5, 5))])
        >>> smoothness_p(10431, power=1)
        (-1, [(3, (2, 2, 2)), (61, (1, 5, 5)), (19, (1, 3, 9))])

    If visual=True then an annotated string will be returned:

        >>> print(smoothness_p(21477639576571, visual=1))
        p**i=4410317**1 has p-1 B=1787, B-pow=1787
        p**i=4869863**1 has p-1 B=2434931, B-pow=2434931

    This string can also be generated directly from a factorization dictionary
    and vice versa:

        >>> factorint(17*9)
        {3: 2, 17: 1}
        >>> smoothness_p(_)
        'p**i=3**2 has p-1 B=2, B-pow=2\\np**i=17**1 has p-1 B=2, B-pow=16'
        >>> smoothness_p(_)
        {3: 2, 17: 1}

    The table of the output logic is:

        ====== ====== ======= =======
        |              Visual
        ------ ----------------------
        Input  True   False   other
        ====== ====== ======= =======
        dict    str    tuple   str
        str     str    tuple   dict
        tuple   str    tuple   str
        n       str    tuple   tuple
        mul     str    tuple   tuple
        ====== ====== ======= =======

    See Also
    ========

    factorint, smoothness
    """

    # visual must be True, False or other (stored as None)
    if visual in (1, 0):
        visual = bool(visual)
    elif visual not in (True, False):
        visual = None

    if isinstance(n, str):
        if visual:
            return n
        d = {}
        for li in n.splitlines():
            k, v = [int(i) for i in
                    li.split('has')[0].split('=')[1].split('**')]
            d[k] = v
        if visual is not True and visual is not False:
            return d
        return smoothness_p(d, visual=False)
    elif not isinstance(n, tuple):
        facs = factorint(n, visual=False)

    if power:
        k = -1
    else:
        k = 1
    if isinstance(n, tuple):
        rv = n
    else:
        rv = (m, sorted([(f,
                          tuple([M] + list(smoothness(f + m))))
                         for f, M in list(facs.items())],
                        key=lambda x: (x[1][k], x[0])))

    if visual is False or (visual is not True) and (type(n) in [int, Mul]):
        return rv
    lines = []
    for dat in rv[1]:
        dat = flatten(dat)
        dat.insert(2, m)
        lines.append('p**i=%i**%i has p%+i B=%i, B-pow=%i' % tuple(dat))
    return '\n'.join(lines)


def trailing(n):
    """Count the number of trailing zero digits in the binary
    representation of n, i.e. determine the largest power of 2
    that divides n.

    Examples
    ========

    >>> from sympy import trailing
    >>> trailing(128)
    7
    >>> trailing(63)
    0
    """
    n = abs(int(n))
    if not n:
        return 0
    low_byte = n & 0xff
    if low_byte:
        return small_trailing[low_byte]

    # 2**m is quick for z up through 2**30
    z = bitcount(n) - 1
    if isinstance(z, SYMPY_INTS):
        if n == 1 << z:
            return z

    if z < 300:
        # fixed 8-byte reduction
        t = 8
        n >>= 8
        while not n & 0xff:
            n >>= 8
            t += 8
        return t + small_trailing[n & 0xff]

    # binary reduction important when there might be a large
    # number of trailing 0s
    t = 0
    p = 8
    while not n & 1:
        while not n & ((1 << p) - 1):
            n >>= p
            t += p
            p *= 2
        p //= 2
    return t


def multiplicity(p, n):
    """
    Find the greatest integer m such that p**m divides n.

    Examples
    ========

    >>> from sympy import multiplicity, Rational
    >>> [multiplicity(5, n) for n in [8, 5, 25, 125, 250]]
    [0, 1, 2, 3, 3]
    >>> multiplicity(3, Rational(1, 9))
    -2

    Note: when checking for the multiplicity of a number in a
    large factorial it is most efficient to send it as an unevaluated
    factorial or to call ``multiplicity_in_factorial`` directly:

    >>> from sympy.ntheory import multiplicity_in_factorial
    >>> from sympy import factorial
    >>> p = factorial(25)
    >>> n = 2**100
    >>> nfac = factorial(n, evaluate=False)
    >>> multiplicity(p, nfac)
    52818775009509558395695966887
    >>> _ == multiplicity_in_factorial(p, n)
    True

    """
    try:
        p, n = as_int(p), as_int(n)
    except ValueError:
        from sympy.functions.combinatorial.factorials import factorial
        if all(isinstance(i, (SYMPY_INTS, Rational)) for i in (p, n)):
            p = Rational(p)
            n = Rational(n)
            if p.q == 1:
                if n.p == 1:
                    return -multiplicity(p.p, n.q)
                return multiplicity(p.p, n.p) - multiplicity(p.p, n.q)
            elif p.p == 1:
                return multiplicity(p.q, n.q)
            else:
                like = min(
                    multiplicity(p.p, n.p),
                    multiplicity(p.q, n.q))
                cross = min(
                    multiplicity(p.q, n.p),
                    multiplicity(p.p, n.q))
                return like - cross
        elif (isinstance(p, (SYMPY_INTS, Integer)) and
                isinstance(n, factorial) and
                isinstance(n.args[0], Integer) and
                n.args[0] >= 0):
            return multiplicity_in_factorial(p, n.args[0])
        raise ValueError('expecting ints or fractions, got %s and %s' % (p, n))

    if n == 0:
        raise ValueError('no such integer exists: multiplicity of %s is not-defined' %(n))
    if p == 2:
        return trailing(n)
    if p < 2:
        raise ValueError('p must be an integer, 2 or larger, but got %s' % p)
    if p == n:
        return 1

    m = 0
    n, rem = divmod(n, p)
    while not rem:
        m += 1
        if m > 5:
            # The multiplicity could be very large. Better
            # to increment in powers of two
            e = 2
            while 1:
                ppow = p**e
                if ppow < n:
                    nnew, rem = divmod(n, ppow)
                    if not rem:
                        m += e
                        e *= 2
                        n = nnew
                        continue
                return m + multiplicity(p, n)
        n, rem = divmod(n, p)
    return m


def multiplicity_in_factorial(p, n):
    """return the largest integer ``m`` such that ``p**m`` divides ``n!``
    without calculating the factorial of ``n``.


    Examples
    ========

    >>> from sympy.ntheory import multiplicity_in_factorial
    >>> from sympy import factorial

    >>> multiplicity_in_factorial(2, 3)
    1

    An instructive use of this is to tell how many trailing zeros
    a given factorial has. For example, there are 6 in 25!:

    >>> factorial(25)
    15511210043330985984000000
    >>> multiplicity_in_factorial(10, 25)
    6

    For large factorials, it is much faster/feasible to use
    this function rather than computing the actual factorial:

    >>> multiplicity_in_factorial(factorial(25), 2**100)
    52818775009509558395695966887

    """

    p, n = as_int(p), as_int(n)

    if p <= 0:
        raise ValueError('expecting positive integer got %s' % p )

    if n < 0:
        raise ValueError('expecting non-negative integer got %s' % n )

    factors = factorint(p)

    # keep only the largest of a given multiplicity since those
    # of a given multiplicity will be goverened by the behavior
    # of the largest factor
    test = defaultdict(int)
    for k, v in factors.items():
        test[v] = max(k, test[v])
    keep = set(test.values())
    # remove others from factors
    for k in list(factors.keys()):
        if k not in keep:
            factors.pop(k)

    mp = S.Infinity
    for i in factors:
        # multiplicity of i in n! is
        mi = (n - (sum(digits(n, i)) - i))//(i - 1)
        # multiplicity of p in n! depends on multiplicity
        # of prime `i` in p, so we floor divide by factors[i]
        # and keep it if smaller than the multiplicity of p
        # seen so far
        mp = min(mp, mi//factors[i])

    return mp


def perfect_power(n, candidates=None, big=True, factor=True):
    """
    Return ``(b, e)`` such that ``n`` == ``b**e`` if ``n`` is a unique
    perfect power with ``e > 1``, else ``False`` (e.g. 1 is not a
    perfect power). A ValueError is raised if ``n`` is not Rational.

    By default, the base is recursively decomposed and the exponents
    collected so the largest possible ``e`` is sought. If ``big=False``
    then the smallest possible ``e`` (thus prime) will be chosen.

    If ``factor=True`` then simultaneous factorization of ``n`` is
    attempted since finding a factor indicates the only possible root
    for ``n``. This is True by default since only a few small factors will
    be tested in the course of searching for the perfect power.

    The use of ``candidates`` is primarily for internal use; if provided,
    False will be returned if ``n`` cannot be written as a power with one
    of the candidates as an exponent and factoring (beyond testing for
    a factor of 2) will not be attempted.

    Examples
    ========

    >>> from sympy import perfect_power, Rational
    >>> perfect_power(16)
    (2, 4)
    >>> perfect_power(16, big=False)
    (4, 2)

    Negative numbers can only have odd perfect powers:

    >>> perfect_power(-4)
    False
    >>> perfect_power(-8)
    (-2, 3)

    Rationals are also recognized:

    >>> perfect_power(Rational(1, 2)**3)
    (1/2, 3)
    >>> perfect_power(Rational(-3, 2)**3)
    (-3/2, 3)

    Notes
    =====

    To know whether an integer is a perfect power of 2 use

        >>> is2pow = lambda n: bool(n and not n & (n - 1))
        >>> [(i, is2pow(i)) for i in range(5)]
        [(0, False), (1, True), (2, True), (3, False), (4, True)]

    It is not necessary to provide ``candidates``. When provided
    it will be assumed that they are ints. The first one that is
    larger than the computed maximum possible exponent will signal
    failure for the routine.

        >>> perfect_power(3**8, [9])
        False
        >>> perfect_power(3**8, [2, 4, 8])
        (3, 8)
        >>> perfect_power(3**8, [4, 8], big=False)
        (9, 4)

    See Also
    ========
    sympy.core.power.integer_nthroot
    sympy.ntheory.primetest.is_square
    """
    if isinstance(n, Rational) and not n.is_Integer:
        p, q = n.as_numer_denom()
        if p is S.One:
            pp = perfect_power(q)
            if pp:
                pp = (n.func(1, pp[0]), pp[1])
        else:
            pp = perfect_power(p)
            if pp:
                num, e = pp
                pq = perfect_power(q, [e])
                if pq:
                    den, _ = pq
                    pp = n.func(num, den), e
        return pp

    n = as_int(n)
    if n < 0:
        pp = perfect_power(-n)
        if pp:
            b, e = pp
            if e % 2:
                return -b, e
        return False

    if n <= 3:
        # no unique exponent for 0, 1
        # 2 and 3 have exponents of 1
        return False
    logn = math.log(n, 2)
    max_possible = int(logn) + 2  # only check values less than this
    not_square = n % 10 in [2, 3, 7, 8]  # squares cannot end in 2, 3, 7, 8
    min_possible = 2 + not_square
    if not candidates:
        candidates = primerange(min_possible, max_possible)
    else:
        candidates = sorted([i for i in candidates
            if min_possible <= i < max_possible])
        if n%2 == 0:
            e = trailing(n)
            candidates = [i for i in candidates if e%i == 0]
        if big:
            candidates = reversed(candidates)
        for e in candidates:
            r, ok = integer_nthroot(n, e)
            if ok:
                return (r, e)
        return False

    def _factors():
        rv = 2 + n % 2
        while True:
            yield rv
            rv = nextprime(rv)

    for fac, e in zip(_factors(), candidates):
        # see if there is a factor present
        if factor and n % fac == 0:
            # find what the potential power is
            if fac == 2:
                e = trailing(n)
            else:
                e = multiplicity(fac, n)
            # if it's a trivial power we are done
            if e == 1:
                return False

            # maybe the e-th root of n is exact
            r, exact = integer_nthroot(n, e)
            if not exact:
                # Having a factor, we know that e is the maximal
                # possible value for a root of n.
                # If n = fac**e*m can be written as a perfect
                # power then see if m can be written as r**E where
                # gcd(e, E) != 1 so n = (fac**(e//E)*r)**E
                m = n//fac**e
                rE = perfect_power(m, candidates=divisors(e, generator=True))
                if not rE:
                    return False
                else:
                    r, E = rE
                    r, e = fac**(e//E)*r, E
            if not big:
                e0 = primefactors(e)
                if e0[0] != e:
                    r, e = r**(e//e0[0]), e0[0]
            return r, e

        # Weed out downright impossible candidates
        if logn/e < 40:
            b = 2.0**(logn/e)
            if abs(int(b + 0.5) - b) > 0.01:
                continue

        # now see if the plausible e makes a perfect power
        r, exact = integer_nthroot(n, e)
        if exact:
            if big:
                m = perfect_power(r, big=big, factor=factor)
                if m:
                    r, e = m[0], e*m[1]
            return int(r), e

    return False


def pollard_rho(n, s=2, a=1, retries=5, seed=1234, max_steps=None, F=None):
    r"""
    Use Pollard's rho method to try to extract a nontrivial factor
    of ``n``. The returned factor may be a composite number. If no
    factor is found, ``None`` is returned.

    The algorithm generates pseudo-random values of x with a generator
    function, replacing x with F(x). If F is not supplied then the
    function x**2 + ``a`` is used. The first value supplied to F(x) is ``s``.
    Upon failure (if ``retries`` is > 0) a new ``a`` and ``s`` will be
    supplied; the ``a`` will be ignored if F was supplied.

    The sequence of numbers generated by such functions generally have a
    a lead-up to some number and then loop around back to that number and
    begin to repeat the sequence, e.g. 1, 2, 3, 4, 5, 3, 4, 5 -- this leader
    and loop look a bit like the Greek letter rho, and thus the name, 'rho'.

    For a given function, very different leader-loop values can be obtained
    so it is a good idea to allow for retries:

    >>> from sympy.ntheory.generate import cycle_length
    >>> n = 16843009
    >>> F = lambda x:(2048*pow(x, 2, n) + 32767) % n
    >>> for s in range(5):
    ...     print('loop length = %4i; leader length = %3i' % next(cycle_length(F, s)))
    ...
    loop length = 2489; leader length =  42
    loop length =   78; leader length = 120
    loop length = 1482; leader length =  99
    loop length = 1482; leader length = 285
    loop length = 1482; leader length = 100

    Here is an explicit example where there is a two element leadup to
    a sequence of 3 numbers (11, 14, 4) that then repeat:

    >>> x=2
    >>> for i in range(9):
    ...     x=(x**2+12)%17
    ...     print(x)
    ...
    16
    13
    11
    14
    4
    11
    14
    4
    11
    >>> next(cycle_length(lambda x: (x**2+12)%17, 2))
    (3, 2)
    >>> list(cycle_length(lambda x: (x**2+12)%17, 2, values=True))
    [16, 13, 11, 14, 4]

    Instead of checking the differences of all generated values for a gcd
    with n, only the kth and 2*kth numbers are checked, e.g. 1st and 2nd,
    2nd and 4th, 3rd and 6th until it has been detected that the loop has been
    traversed. Loops may be many thousands of steps long before rho finds a
    factor or reports failure. If ``max_steps`` is specified, the iteration
    is cancelled with a failure after the specified number of steps.

    Examples
    ========

    >>> from sympy import pollard_rho
    >>> n=16843009
    >>> F=lambda x:(2048*pow(x,2,n) + 32767) % n
    >>> pollard_rho(n, F=F)
    257

    Use the default setting with a bad value of ``a`` and no retries:

    >>> pollard_rho(n, a=n-2, retries=0)

    If retries is > 0 then perhaps the problem will correct itself when
    new values are generated for a:

    >>> pollard_rho(n, a=n-2, retries=1)
    257

    References
    ==========

    .. [1] Richard Crandall & Carl Pomerance (2005), "Prime Numbers:
           A Computational Perspective", Springer, 2nd edition, 229-231

    """
    n = int(n)
    if n < 5:
        raise ValueError('pollard_rho should receive n > 4')
    prng = random.Random(seed + retries)
    V = s
    for i in range(retries + 1):
        U = V
        if not F:
            F = lambda x: (pow(x, 2, n) + a) % n
        j = 0
        while 1:
            if max_steps and (j > max_steps):
                break
            j += 1
            U = F(U)
            V = F(F(V))  # V is 2x further along than U
            g = igcd(U - V, n)
            if g == 1:
                continue
            if g == n:
                break
            return int(g)
        V = prng.randint(0, n - 1)
        a = prng.randint(1, n - 3)  # for x**2 + a, a%n should not be 0 or -2
        F = None
    return None


def pollard_pm1(n, B=10, a=2, retries=0, seed=1234):
    """
    Use Pollard's p-1 method to try to extract a nontrivial factor
    of ``n``. Either a divisor (perhaps composite) or ``None`` is returned.

    The value of ``a`` is the base that is used in the test gcd(a**M - 1, n).
    The default is 2.  If ``retries`` > 0 then if no factor is found after the
    first attempt, a new ``a`` will be generated randomly (using the ``seed``)
    and the process repeated.

    Note: the value of M is lcm(1..B) = reduce(ilcm, range(2, B + 1)).

    A search is made for factors next to even numbers having a power smoothness
    less than ``B``. Choosing a larger B increases the likelihood of finding a
    larger factor but takes longer. Whether a factor of n is found or not
    depends on ``a`` and the power smoothness of the even number just less than
    the factor p (hence the name p - 1).

    Although some discussion of what constitutes a good ``a`` some
    descriptions are hard to interpret. At the modular.math site referenced
    below it is stated that if gcd(a**M - 1, n) = N then a**M % q**r is 1
    for every prime power divisor of N. But consider the following:

        >>> from sympy.ntheory.factor_ import smoothness_p, pollard_pm1
        >>> n=257*1009
        >>> smoothness_p(n)
        (-1, [(257, (1, 2, 256)), (1009, (1, 7, 16))])

    So we should (and can) find a root with B=16:

        >>> pollard_pm1(n, B=16, a=3)
        1009

    If we attempt to increase B to 256 we find that it does not work:

        >>> pollard_pm1(n, B=256)
        >>>

    But if the value of ``a`` is changed we find that only multiples of
    257 work, e.g.:

        >>> pollard_pm1(n, B=256, a=257)
        1009

    Checking different ``a`` values shows that all the ones that did not
    work had a gcd value not equal to ``n`` but equal to one of the
    factors:

        >>> from sympy import ilcm, igcd, factorint, Pow
        >>> M = 1
        >>> for i in range(2, 256):
        ...     M = ilcm(M, i)
        ...
        >>> set([igcd(pow(a, M, n) - 1, n) for a in range(2, 256) if
        ...      igcd(pow(a, M, n) - 1, n) != n])
        {1009}

    But does aM % d for every divisor of n give 1?

        >>> aM = pow(255, M, n)
        >>> [(d, aM%Pow(*d.args)) for d in factorint(n, visual=True).args]
        [(257**1, 1), (1009**1, 1)]

    No, only one of them. So perhaps the principle is that a root will
    be found for a given value of B provided that:

    1) the power smoothness of the p - 1 value next to the root
       does not exceed B
    2) a**M % p != 1 for any of the divisors of n.

    By trying more than one ``a`` it is possible that one of them
    will yield a factor.

    Examples
    ========

    With the default smoothness bound, this number cannot be cracked:

        >>> from sympy.ntheory import pollard_pm1
        >>> pollard_pm1(21477639576571)

    Increasing the smoothness bound helps:

        >>> pollard_pm1(21477639576571, B=2000)
        4410317

    Looking at the smoothness of the factors of this number we find:

        >>> from sympy.ntheory.factor_ import smoothness_p, factorint
        >>> print(smoothness_p(21477639576571, visual=1))
        p**i=4410317**1 has p-1 B=1787, B-pow=1787
        p**i=4869863**1 has p-1 B=2434931, B-pow=2434931

    The B and B-pow are the same for the p - 1 factorizations of the divisors
    because those factorizations had a very large prime factor:

        >>> factorint(4410317 - 1)
        {2: 2, 617: 1, 1787: 1}
        >>> factorint(4869863-1)
        {2: 1, 2434931: 1}

    Note that until B reaches the B-pow value of 1787, the number is not cracked;

        >>> pollard_pm1(21477639576571, B=1786)
        >>> pollard_pm1(21477639576571, B=1787)
        4410317

    The B value has to do with the factors of the number next to the divisor,
    not the divisors themselves. A worst case scenario is that the number next
    to the factor p has a large prime divisisor or is a perfect power. If these
    conditions apply then the power-smoothness will be about p/2 or p. The more
    realistic is that there will be a large prime factor next to p requiring
    a B value on the order of p/2. Although primes may have been searched for
    up to this level, the p/2 is a factor of p - 1, something that we do not
    know. The modular.math reference below states that 15% of numbers in the
    range of 10**15 to 15**15 + 10**4 are 10**6 power smooth so a B of 10**6
    will fail 85% of the time in that range. From 10**8 to 10**8 + 10**3 the
    percentages are nearly reversed...but in that range the simple trial
    division is quite fast.

    References
    ==========

    .. [1] Richard Crandall & Carl Pomerance (2005), "Prime Numbers:
           A Computational Perspective", Springer, 2nd edition, 236-238
    .. [2] https://web.archive.org/web/20150716201437/http://modular.math.washington.edu/edu/2007/spring/ent/ent-html/node81.html
    .. [3] https://www.cs.toronto.edu/~yuvalf/Factorization.pdf
    """

    n = int(n)
    if n < 4 or B < 3:
        raise ValueError('pollard_pm1 should receive n > 3 and B > 2')
    prng = random.Random(seed + B)

    # computing a**lcm(1,2,3,..B) % n for B > 2
    # it looks weird, but it's right: primes run [2, B]
    # and the answer's not right until the loop is done.
    for i in range(retries + 1):
        aM = a
        for p in sieve.primerange(2, B + 1):
            e = int(math.log(B, p))
            aM = pow(aM, pow(p, e), n)
        g = igcd(aM - 1, n)
        if 1 < g < n:
            return int(g)

        # get a new a:
        # since the exponent, lcm(1..B), is even, if we allow 'a' to be 'n-1'
        # then (n - 1)**even % n will be 1 which will give a g of 0 and 1 will
        # give a zero, too, so we set the range as [2, n-2]. Some references
        # say 'a' should be coprime to n, but either will detect factors.
        a = prng.randint(2, n - 2)


def _trial(factors, n, candidates, verbose=False):
    """
    Helper function for integer factorization. Trial factors ``n`
    against all integers given in the sequence ``candidates``
    and updates the dict ``factors`` in-place. Returns the reduced
    value of ``n`` and a flag indicating whether any factors were found.
    """
    if verbose:
        factors0 = list(factors.keys())
    nfactors = len(factors)
    for d in candidates:
        if n % d == 0:
            m = multiplicity(d, n)
            n //= d**m
            factors[d] = m
    if verbose:
        for k in sorted(set(factors).difference(set(factors0))):
            print(factor_msg % (k, factors[k]))
    return int(n), len(factors) != nfactors


def _check_termination(factors, n, limitp1, use_trial, use_rho, use_pm1,
                       verbose):
    """
    Helper function for integer factorization. Checks if ``n``
    is a prime or a perfect power, and in those cases updates
    the factorization and raises ``StopIteration``.
    """

    if verbose:
        print('Check for termination')

    # since we've already been factoring there is no need to do
    # simultaneous factoring with the power check
    p = perfect_power(n, factor=False)
    if p is not False:
        base, exp = p
        if limitp1:
            limit = limitp1 - 1
        else:
            limit = limitp1
        facs = factorint(base, limit, use_trial, use_rho, use_pm1,
                         verbose=False)
        for b, e in facs.items():
            if verbose:
                print(factor_msg % (b, e))
            factors[b] = exp*e
        raise StopIteration

    if isprime(n):
        factors[int(n)] = 1
        raise StopIteration

    if n == 1:
        raise StopIteration

trial_int_msg = "Trial division with ints [%i ... %i] and fail_max=%i"
trial_msg = "Trial division with primes [%i ... %i]"
rho_msg = "Pollard's rho with retries %i, max_steps %i and seed %i"
pm1_msg = "Pollard's p-1 with smoothness bound %i and seed %i"
ecm_msg = "Elliptic Curve with B1 bound %i, B2 bound %i, num_curves %i"
factor_msg = '\t%i ** %i'
fermat_msg = 'Close factors satisying Fermat condition found.'
complete_msg = 'Factorization is complete.'


def _factorint_small(factors, n, limit, fail_max):
    """
    Return the value of n and either a 0 (indicating that factorization up
    to the limit was complete) or else the next near-prime that would have
    been tested.

    Factoring stops if there are fail_max unsuccessful tests in a row.

    If factors of n were found they will be in the factors dictionary as
    {factor: multiplicity} and the returned value of n will have had those
    factors removed. The factors dictionary is modified in-place.

    """

    def done(n, d):
        """return n, d if the sqrt(n) was not reached yet, else
           n, 0 indicating that factoring is done.
        """
        if d*d <= n:
            return n, d
        return n, 0

    d = 2
    m = trailing(n)
    if m:
        factors[d] = m
        n >>= m
    d = 3
    if limit < d:
        if n > 1:
            factors[n] = 1
        return done(n, d)
    # reduce
    m = 0
    while n % d == 0:
        n //= d
        m += 1
        if m == 20:
            mm = multiplicity(d, n)
            m += mm
            n //= d**mm
            break
    if m:
        factors[d] = m

    # when d*d exceeds maxx or n we are done; if limit**2 is greater
    # than n then maxx is set to zero so the value of n will flag the finish
    if limit*limit > n:
        maxx = 0
    else:
        maxx = limit*limit

    dd = maxx or n
    d = 5
    fails = 0
    while fails < fail_max:
        if d*d > dd:
            break
        # d = 6*i - 1
        # reduce
        m = 0
        while n % d == 0:
            n //= d
            m += 1
            if m == 20:
                mm = multiplicity(d, n)
                m += mm
                n //= d**mm
                break
        if m:
            factors[d] = m
            dd = maxx or n
            fails = 0
        else:
            fails += 1
        d += 2
        if d*d > dd:
            break
        # d = 6*i - 1
        # reduce
        m = 0
        while n % d == 0:
            n //= d
            m += 1
            if m == 20:
                mm = multiplicity(d, n)
                m += mm
                n //= d**mm
                break
        if m:
            factors[d] = m
            dd = maxx or n
            fails = 0
        else:
            fails += 1
        # d = 6*(i + 1) - 1
        d += 4

    return done(n, d)


def factorint(n, limit=None, use_trial=True, use_rho=True, use_pm1=True,
              use_ecm=True, verbose=False, visual=None, multiple=False):
    r"""
    Given a positive integer ``n``, ``factorint(n)`` returns a dict containing
    the prime factors of ``n`` as keys and their respective multiplicities
    as values. For example:

    >>> from sympy.ntheory import factorint
    >>> factorint(2000)    # 2000 = (2**4) * (5**3)
    {2: 4, 5: 3}
    >>> factorint(65537)   # This number is prime
    {65537: 1}

    For input less than 2, factorint behaves as follows:

        - ``factorint(1)`` returns the empty factorization, ``{}``
        - ``factorint(0)`` returns ``{0:1}``
        - ``factorint(-n)`` adds ``-1:1`` to the factors and then factors ``n``

    Partial Factorization:

    If ``limit`` (> 3) is specified, the search is stopped after performing
    trial division up to (and including) the limit (or taking a
    corresponding number of rho/p-1 steps). This is useful if one has
    a large number and only is interested in finding small factors (if
    any). Note that setting a limit does not prevent larger factors
    from being found early; it simply means that the largest factor may
    be composite. Since checking for perfect power is relatively cheap, it is
    done regardless of the limit setting.

    This number, for example, has two small factors and a huge
    semi-prime factor that cannot be reduced easily:

    >>> from sympy.ntheory import isprime
    >>> a = 1407633717262338957430697921446883
    >>> f = factorint(a, limit=10000)
    >>> f == {991: 1, int(202916782076162456022877024859): 1, 7: 1}
    True
    >>> isprime(max(f))
    False

    This number has a small factor and a residual perfect power whose
    base is greater than the limit:

    >>> factorint(3*101**7, limit=5)
    {3: 1, 101: 7}

    List of Factors:

    If ``multiple`` is set to ``True`` then a list containing the
    prime factors including multiplicities is returned.

    >>> factorint(24, multiple=True)
    [2, 2, 2, 3]

    Visual Factorization:

    If ``visual`` is set to ``True``, then it will return a visual
    factorization of the integer.  For example:

    >>> from sympy import pprint
    >>> pprint(factorint(4200, visual=True))
     3  1  2  1
    2 *3 *5 *7

    Note that this is achieved by using the evaluate=False flag in Mul
    and Pow. If you do other manipulations with an expression where
    evaluate=False, it may evaluate.  Therefore, you should use the
    visual option only for visualization, and use the normal dictionary
    returned by visual=False if you want to perform operations on the
    factors.

    You can easily switch between the two forms by sending them back to
    factorint:

    >>> from sympy import Mul
    >>> regular = factorint(1764); regular
    {2: 2, 3: 2, 7: 2}
    >>> pprint(factorint(regular))
     2  2  2
    2 *3 *7

    >>> visual = factorint(1764, visual=True); pprint(visual)
     2  2  2
    2 *3 *7
    >>> print(factorint(visual))
    {2: 2, 3: 2, 7: 2}

    If you want to send a number to be factored in a partially factored form
    you can do so with a dictionary or unevaluated expression:

    >>> factorint(factorint({4: 2, 12: 3})) # twice to toggle to dict form
    {2: 10, 3: 3}
    >>> factorint(Mul(4, 12, evaluate=False))
    {2: 4, 3: 1}

    The table of the output logic is:

        ====== ====== ======= =======
                       Visual
        ------ ----------------------
        Input  True   False   other
        ====== ====== ======= =======
        dict    mul    dict    mul
        n       mul    dict    dict
        mul     mul    dict    dict
        ====== ====== ======= =======

    Notes
    =====

    Algorithm:

    The function switches between multiple algorithms. Trial division
    quickly finds small factors (of the order 1-5 digits), and finds
    all large factors if given enough time. The Pollard rho and p-1
    algorithms are used to find large factors ahead of time; they
    will often find factors of the order of 10 digits within a few
    seconds:

    >>> factors = factorint(12345678910111213141516)
    >>> for base, exp in sorted(factors.items()):
    ...     print('%s %s' % (base, exp))
    ...
    2 2
    2507191691 1
    1231026625769 1

    Any of these methods can optionally be disabled with the following
    boolean parameters:

        - ``use_trial``: Toggle use of trial division
        - ``use_rho``: Toggle use of Pollard's rho method
        - ``use_pm1``: Toggle use of Pollard's p-1 method

    ``factorint`` also periodically checks if the remaining part is
    a prime number or a perfect power, and in those cases stops.

    For unevaluated factorial, it uses Legendre's formula(theorem).


    If ``verbose`` is set to ``True``, detailed progress is printed.

    See Also
    ========

    smoothness, smoothness_p, divisors

    """
    if isinstance(n, Dict):
        n = dict(n)
    if multiple:
        fac = factorint(n, limit=limit, use_trial=use_trial,
                           use_rho=use_rho, use_pm1=use_pm1,
                           verbose=verbose, visual=False, multiple=False)
        factorlist = sum(([p] * fac[p] if fac[p] > 0 else [S.One/p]*(-fac[p])
                               for p in sorted(fac)), [])
        return factorlist

    factordict = {}
    if visual and not isinstance(n, (Mul, dict)):
        factordict = factorint(n, limit=limit, use_trial=use_trial,
                               use_rho=use_rho, use_pm1=use_pm1,
                               verbose=verbose, visual=False)
    elif isinstance(n, Mul):
        factordict = {int(k): int(v) for k, v in
            n.as_powers_dict().items()}
    elif isinstance(n, dict):
        factordict = n
    if factordict and isinstance(n, (Mul, dict)):
        # check it
        for key in list(factordict.keys()):
            if isprime(key):
                continue
            e = factordict.pop(key)
            d = factorint(key, limit=limit, use_trial=use_trial, use_rho=use_rho,
                          use_pm1=use_pm1, verbose=verbose, visual=False)
            for k, v in d.items():
                if k in factordict:
                    factordict[k] += v*e
                else:
                    factordict[k] = v*e
    if visual or (type(n) is dict and
                  visual is not True and
                  visual is not False):
        if factordict == {}:
            return S.One
        if -1 in factordict:
            factordict.pop(-1)
            args = [S.NegativeOne]
        else:
            args = []
        args.extend([Pow(*i, evaluate=False)
                     for i in sorted(factordict.items())])
        return Mul(*args, evaluate=False)
    elif isinstance(n, (dict, Mul)):
        return factordict

    assert use_trial or use_rho or use_pm1 or use_ecm

    from sympy.functions.combinatorial.factorials import factorial
    if isinstance(n, factorial):
        x = as_int(n.args[0])
        if x >= 20:
            factors = {}
            m = 2 # to initialize the if condition below
            for p in sieve.primerange(2, x + 1):
                if m > 1:
                    m, q = 0, x // p
                    while q != 0:
                        m += q
                        q //= p
                factors[p] = m
            if factors and verbose:
                for k in sorted(factors):
                    print(factor_msg % (k, factors[k]))
            if verbose:
                print(complete_msg)
            return factors
        else:
            # if n < 20!, direct computation is faster
            # since it uses a lookup table
            n = n.func(x)

    n = as_int(n)
    if limit:
        limit = int(limit)
        use_ecm = False

    # special cases
    if n < 0:
        factors = factorint(
            -n, limit=limit, use_trial=use_trial, use_rho=use_rho,
            use_pm1=use_pm1, verbose=verbose, visual=False)
        factors[-1] = 1
        return factors

    if limit and limit < 2:
        if n == 1:
            return {}
        return {n: 1}
    elif n < 10:
        # doing this we are assured of getting a limit > 2
        # when we have to compute it later
        return [{0: 1}, {}, {2: 1}, {3: 1}, {2: 2}, {5: 1},
                {2: 1, 3: 1}, {7: 1}, {2: 3}, {3: 2}][n]

    factors = {}

    # do simplistic factorization
    if verbose:
        sn = str(n)
        if len(sn) > 50:
            print('Factoring %s' % sn[:5] + \
                  '..(%i other digits)..' % (len(sn) - 10) + sn[-5:])
        else:
            print('Factoring', n)

    if use_trial:
        # this is the preliminary factorization for small factors
        small = 2**15
        fail_max = 600
        small = min(small, limit or small)
        if verbose:
            print(trial_int_msg % (2, small, fail_max))
        n, next_p = _factorint_small(factors, n, small, fail_max)
    else:
        next_p = 2
    if factors and verbose:
        for k in sorted(factors):
            print(factor_msg % (k, factors[k]))
    if next_p == 0:
        if n > 1:
            factors[int(n)] = 1
        if verbose:
            print(complete_msg)
        return factors

    # continue with more advanced factorization methods

    # first check if the simplistic run didn't finish
    # because of the limit and check for a perfect
    # power before exiting
    try:
        if limit and next_p > limit:
            if verbose:
                print('Exceeded limit:', limit)

            _check_termination(factors, n, limit, use_trial, use_rho, use_pm1,
                               verbose)

            if n > 1:
                factors[int(n)] = 1
            return factors
        else:
            # Before quitting (or continuing on)...

            # ...do a Fermat test since it's so easy and we need the
            # square root anyway. Finding 2 factors is easy if they are
            # "close enough." This is the big root equivalent of dividing by
            # 2, 3, 5.
            sqrt_n = integer_nthroot(n, 2)[0]
            a = sqrt_n + 1
            a2 = a**2
            b2 = a2 - n
            for i in range(3):
                b, fermat = integer_nthroot(b2, 2)
                if fermat:
                    break
                b2 += 2*a + 1  # equiv to (a + 1)**2 - n
                a += 1
            if fermat:
                if verbose:
                    print(fermat_msg)
                if limit:
                    limit -= 1
                for r in [a - b, a + b]:
                    facs = factorint(r, limit=limit, use_trial=use_trial,
                                     use_rho=use_rho, use_pm1=use_pm1,
                                     verbose=verbose)
                    for k, v in facs.items():
                        factors[k] = factors.get(k, 0) + v
                raise StopIteration

            # ...see if factorization can be terminated
            _check_termination(factors, n, limit, use_trial, use_rho, use_pm1,
                               verbose)

    except StopIteration:
        if verbose:
            print(complete_msg)
        return factors

    # these are the limits for trial division which will
    # be attempted in parallel with pollard methods
    low, high = next_p, 2*next_p

    limit = limit or sqrt_n
    # add 1 to make sure limit is reached in primerange calls
    limit += 1
    iteration = 0
    while 1:

        try:
            high_ = high
            if limit < high_:
                high_ = limit

            # Trial division
            if use_trial:
                if verbose:
                    print(trial_msg % (low, high_))
                ps = sieve.primerange(low, high_)
                n, found_trial = _trial(factors, n, ps, verbose)
                if found_trial:
                    _check_termination(factors, n, limit, use_trial, use_rho,
                                       use_pm1, verbose)
            else:
                found_trial = False

            if high > limit:
                if verbose:
                    print('Exceeded limit:', limit)
                if n > 1:
                    factors[int(n)] = 1
                raise StopIteration

            # Only used advanced methods when no small factors were found
            if not found_trial:
                if (use_pm1 or use_rho):
                    high_root = max(int(math.log(high_**0.7)), low, 3)

                    # Pollard p-1
                    if use_pm1:
                        if verbose:
                            print(pm1_msg % (high_root, high_))
                        c = pollard_pm1(n, B=high_root, seed=high_)
                        if c:
                            # factor it and let _trial do the update
                            ps = factorint(c, limit=limit - 1,
                                           use_trial=use_trial,
                                           use_rho=use_rho,
                                           use_pm1=use_pm1,
                                           use_ecm=use_ecm,
                                           verbose=verbose)
                            n, _ = _trial(factors, n, ps, verbose=False)
                            _check_termination(factors, n, limit, use_trial,
                                               use_rho, use_pm1, verbose)

                    # Pollard rho
                    if use_rho:
                        max_steps = high_root
                        if verbose:
                            print(rho_msg % (1, max_steps, high_))
                        c = pollard_rho(n, retries=1, max_steps=max_steps,
                                        seed=high_)
                        if c:
                            # factor it and let _trial do the update
                            ps = factorint(c, limit=limit - 1,
                                           use_trial=use_trial,
                                           use_rho=use_rho,
                                           use_pm1=use_pm1,
                                           use_ecm=use_ecm,
                                           verbose=verbose)
                            n, _ = _trial(factors, n, ps, verbose=False)
                            _check_termination(factors, n, limit, use_trial,
                                               use_rho, use_pm1, verbose)

        except StopIteration:
            if verbose:
                print(complete_msg)
            return factors
        #Use subexponential algorithms if use_ecm
        #Use pollard algorithms for finding small factors for 3 iterations
        #if after small factors the number of digits of n is >= 20 then use ecm
        iteration += 1
        if use_ecm and iteration >= 3 and len(str(n)) >= 25:
            break
        low, high = high, high*2
    B1 = 10000
    B2 = 100*B1
    num_curves = 50
    while(1):
        if verbose:
            print(ecm_msg % (B1, B2, num_curves))
        while(1):
            try:
                factor = _ecm_one_factor(n, B1, B2, num_curves)
                ps = factorint(factor, limit=limit - 1,
                               use_trial=use_trial,
                               use_rho=use_rho,
                               use_pm1=use_pm1,
                               use_ecm=use_ecm,
                               verbose=verbose)
                n, _ = _trial(factors, n, ps, verbose=False)
                _check_termination(factors, n, limit, use_trial,
                                       use_rho, use_pm1, verbose)
            except ValueError:
                break
            except StopIteration:
                if verbose:
                    print(complete_msg)
                return factors
        B1 *= 5
        B2 = 100*B1
        num_curves *= 4


def factorrat(rat, limit=None, use_trial=True, use_rho=True, use_pm1=True,
              verbose=False, visual=None, multiple=False):
    r"""
    Given a Rational ``r``, ``factorrat(r)`` returns a dict containing
    the prime factors of ``r`` as keys and their respective multiplicities
    as values. For example:

    >>> from sympy import factorrat, S
    >>> factorrat(S(8)/9)    # 8/9 = (2**3) * (3**-2)
    {2: 3, 3: -2}
    >>> factorrat(S(-1)/987)    # -1/789 = -1 * (3**-1) * (7**-1) * (47**-1)
    {-1: 1, 3: -1, 7: -1, 47: -1}

    Please see the docstring for ``factorint`` for detailed explanations
    and examples of the following keywords:

        - ``limit``: Integer limit up to which trial division is done
        - ``use_trial``: Toggle use of trial division
        - ``use_rho``: Toggle use of Pollard's rho method
        - ``use_pm1``: Toggle use of Pollard's p-1 method
        - ``verbose``: Toggle detailed printing of progress
        - ``multiple``: Toggle returning a list of factors or dict
        - ``visual``: Toggle product form of output
    """
    if multiple:
        fac = factorrat(rat, limit=limit, use_trial=use_trial,
                  use_rho=use_rho, use_pm1=use_pm1,
                  verbose=verbose, visual=False, multiple=False)
        factorlist = sum(([p] * fac[p] if fac[p] > 0 else [S.One/p]*(-fac[p])
                               for p, _ in sorted(fac.items(),
                                                        key=lambda elem: elem[0]
                                                        if elem[1] > 0
                                                        else 1/elem[0])), [])
        return factorlist

    f = factorint(rat.p, limit=limit, use_trial=use_trial,
                  use_rho=use_rho, use_pm1=use_pm1,
                  verbose=verbose).copy()
    f = defaultdict(int, f)
    for p, e in factorint(rat.q, limit=limit,
                          use_trial=use_trial,
                          use_rho=use_rho,
                          use_pm1=use_pm1,
                          verbose=verbose).items():
        f[p] += -e

    if len(f) > 1 and 1 in f:
        del f[1]
    if not visual:
        return dict(f)
    else:
        if -1 in f:
            f.pop(-1)
            args = [S.NegativeOne]
        else:
            args = []
        args.extend([Pow(*i, evaluate=False)
                     for i in sorted(f.items())])
        return Mul(*args, evaluate=False)



def primefactors(n, limit=None, verbose=False):
    """Return a sorted list of n's prime factors, ignoring multiplicity
    and any composite factor that remains if the limit was set too low
    for complete factorization. Unlike factorint(), primefactors() does
    not return -1 or 0.

    Examples
    ========

    >>> from sympy.ntheory import primefactors, factorint, isprime
    >>> primefactors(6)
    [2, 3]
    >>> primefactors(-5)
    [5]

    >>> sorted(factorint(123456).items())
    [(2, 6), (3, 1), (643, 1)]
    >>> primefactors(123456)
    [2, 3, 643]

    >>> sorted(factorint(10000000001, limit=200).items())
    [(101, 1), (99009901, 1)]
    >>> isprime(99009901)
    False
    >>> primefactors(10000000001, limit=300)
    [101]

    See Also
    ========

    divisors
    """
    n = int(n)
    factors = sorted(factorint(n, limit=limit, verbose=verbose).keys())
    s = [f for f in factors[:-1:] if f not in [-1, 0, 1]]
    if factors and isprime(factors[-1]):
        s += [factors[-1]]
    return s


def _divisors(n, proper=False):
    """Helper function for divisors which generates the divisors."""

    factordict = factorint(n)
    ps = sorted(factordict.keys())

    def rec_gen(n=0):
        if n == len(ps):
            yield 1
        else:
            pows = [1]
            for j in range(factordict[ps[n]]):
                pows.append(pows[-1] * ps[n])
            for q in rec_gen(n + 1):
                for p in pows:
                    yield p * q

    if proper:
        for p in rec_gen():
            if p != n:
                yield p
    else:
        yield from rec_gen()


def divisors(n, generator=False, proper=False):
    r"""
    Return all divisors of n sorted from 1..n by default.
    If generator is ``True`` an unordered generator is returned.

    The number of divisors of n can be quite large if there are many
    prime factors (counting repeated factors). If only the number of
    factors is desired use divisor_count(n).

    Examples
    ========

    >>> from sympy import divisors, divisor_count
    >>> divisors(24)
    [1, 2, 3, 4, 6, 8, 12, 24]
    >>> divisor_count(24)
    8

    >>> list(divisors(120, generator=True))
    [1, 2, 4, 8, 3, 6, 12, 24, 5, 10, 20, 40, 15, 30, 60, 120]

    Notes
    =====

    This is a slightly modified version of Tim Peters referenced at:
    https://stackoverflow.com/questions/1010381/python-factorization

    See Also
    ========

    primefactors, factorint, divisor_count
    """

    n = as_int(abs(n))
    if isprime(n):
        if proper:
            return [1]
        return [1, n]
    if n == 1:
        if proper:
            return []
        return [1]
    if n == 0:
        return []
    rv = _divisors(n, proper)
    if not generator:
        return sorted(rv)
    return rv


def divisor_count(n, modulus=1, proper=False):
    """
    Return the number of divisors of ``n``. If ``modulus`` is not 1 then only
    those that are divisible by ``modulus`` are counted. If ``proper`` is True
    then the divisor of ``n`` will not be counted.

    Examples
    ========

    >>> from sympy import divisor_count
    >>> divisor_count(6)
    4
    >>> divisor_count(6, 2)
    2
    >>> divisor_count(6, proper=True)
    3

    See Also
    ========

    factorint, divisors, totient, proper_divisor_count

    """

    if not modulus:
        return 0
    elif modulus != 1:
        n, r = divmod(n, modulus)
        if r:
            return 0
    if n == 0:
        return 0
    n = Mul(*[v + 1 for k, v in factorint(n).items() if k > 1])
    if n and proper:
        n -= 1
    return n


def proper_divisors(n, generator=False):
    """
    Return all divisors of n except n, sorted by default.
    If generator is ``True`` an unordered generator is returned.

    Examples
    ========

    >>> from sympy import proper_divisors, proper_divisor_count
    >>> proper_divisors(24)
    [1, 2, 3, 4, 6, 8, 12]
    >>> proper_divisor_count(24)
    7
    >>> list(proper_divisors(120, generator=True))
    [1, 2, 4, 8, 3, 6, 12, 24, 5, 10, 20, 40, 15, 30, 60]

    See Also
    ========

    factorint, divisors, proper_divisor_count

    """
    return divisors(n, generator=generator, proper=True)


def proper_divisor_count(n, modulus=1):
    """
    Return the number of proper divisors of ``n``.

    Examples
    ========

    >>> from sympy import proper_divisor_count
    >>> proper_divisor_count(6)
    3
    >>> proper_divisor_count(6, modulus=2)
    1

    See Also
    ========

    divisors, proper_divisors, divisor_count

    """
    return divisor_count(n, modulus=modulus, proper=True)


def _udivisors(n):
    """Helper function for udivisors which generates the unitary divisors."""

    factorpows = [p**e for p, e in factorint(n).items()]
    for i in range(2**len(factorpows)):
        d, j, k = 1, i, 0
        while j:
            if (j & 1):
                d *= factorpows[k]
            j >>= 1
            k += 1
        yield d


def udivisors(n, generator=False):
    r"""
    Return all unitary divisors of n sorted from 1..n by default.
    If generator is ``True`` an unordered generator is returned.

    The number of unitary divisors of n can be quite large if there are many
    prime factors. If only the number of unitary divisors is desired use
    udivisor_count(n).

    Examples
    ========

    >>> from sympy.ntheory.factor_ import udivisors, udivisor_count
    >>> udivisors(15)
    [1, 3, 5, 15]
    >>> udivisor_count(15)
    4

    >>> sorted(udivisors(120, generator=True))
    [1, 3, 5, 8, 15, 24, 40, 120]

    See Also
    ========

    primefactors, factorint, divisors, divisor_count, udivisor_count

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Unitary_divisor
    .. [2] https://mathworld.wolfram.com/UnitaryDivisor.html

    """

    n = as_int(abs(n))
    if isprime(n):
        return [1, n]
    if n == 1:
        return [1]
    if n == 0:
        return []
    rv = _udivisors(n)
    if not generator:
        return sorted(rv)
    return rv


def udivisor_count(n):
    """
    Return the number of unitary divisors of ``n``.

    Parameters
    ==========

    n : integer

    Examples
    ========

    >>> from sympy.ntheory.factor_ import udivisor_count
    >>> udivisor_count(120)
    8

    See Also
    ========

    factorint, divisors, udivisors, divisor_count, totient

    References
    ==========

    .. [1] https://mathworld.wolfram.com/UnitaryDivisorFunction.html

    """

    if n == 0:
        return 0
    return 2**len([p for p in factorint(n) if p > 1])


def _antidivisors(n):
    """Helper function for antidivisors which generates the antidivisors."""

    for d in _divisors(n):
        y = 2*d
        if n > y and n % y:
            yield y
    for d in _divisors(2*n-1):
        if n > d >= 2 and n % d:
            yield d
    for d in _divisors(2*n+1):
        if n > d >= 2 and n % d:
            yield d


def antidivisors(n, generator=False):
    r"""
    Return all antidivisors of n sorted from 1..n by default.

    Antidivisors [1]_ of n are numbers that do not divide n by the largest
    possible margin.  If generator is True an unordered generator is returned.

    Examples
    ========

    >>> from sympy.ntheory.factor_ import antidivisors
    >>> antidivisors(24)
    [7, 16]

    >>> sorted(antidivisors(128, generator=True))
    [3, 5, 15, 17, 51, 85]

    See Also
    ========

    primefactors, factorint, divisors, divisor_count, antidivisor_count

    References
    ==========

    .. [1] definition is described in https://oeis.org/A066272/a066272a.html

    """

    n = as_int(abs(n))
    if n <= 2:
        return []
    rv = _antidivisors(n)
    if not generator:
        return sorted(rv)
    return rv


def antidivisor_count(n):
    """
    Return the number of antidivisors [1]_ of ``n``.

    Parameters
    ==========

    n : integer

    Examples
    ========

    >>> from sympy.ntheory.factor_ import antidivisor_count
    >>> antidivisor_count(13)
    4
    >>> antidivisor_count(27)
    5

    See Also
    ========

    factorint, divisors, antidivisors, divisor_count, totient

    References
    ==========

    .. [1] formula from https://oeis.org/A066272

    """

    n = as_int(abs(n))
    if n <= 2:
        return 0
    return divisor_count(2*n - 1) + divisor_count(2*n + 1) + \
        divisor_count(n) - divisor_count(n, 2) - 5


class totient(Function):
    r"""
    Calculate the Euler totient function phi(n)

    ``totient(n)`` or `\phi(n)` is the number of positive integers `\leq` n
    that are relatively prime to n.

    Parameters
    ==========

    n : integer

    Examples
    ========

    >>> from sympy.ntheory import totient
    >>> totient(1)
    1
    >>> totient(25)
    20
    >>> totient(45) == totient(5)*totient(9)
    True

    See Also
    ========

    divisor_count

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Euler%27s_totient_function
    .. [2] https://mathworld.wolfram.com/TotientFunction.html

    """
    @classmethod
    def eval(cls, n):
        if n.is_Integer:
            if n < 1:
                raise ValueError("n must be a positive integer")
            factors = factorint(n)
            return cls._from_factors(factors)
        elif not isinstance(n, Expr) or (n.is_integer is False) or (n.is_positive is False):
            raise ValueError("n must be a positive integer")

    def _eval_is_integer(self):
        return fuzzy_and([self.args[0].is_integer, self.args[0].is_positive])

    @classmethod
    def _from_distinct_primes(self, *args):
        """Subroutine to compute totient from the list of assumed
        distinct primes

        Examples
        ========

        >>> from sympy.ntheory.factor_ import totient
        >>> totient._from_distinct_primes(5, 7)
        24
        """
        return reduce(lambda i, j: i * (j-1), args, 1)

    @classmethod
    def _from_factors(self, factors):
        """Subroutine to compute totient from already-computed factors

        Examples
        ========

        >>> from sympy.ntheory.factor_ import totient
        >>> totient._from_factors({5: 2})
        20
        """
        t = 1
        for p, k in factors.items():
            t *= (p - 1) * p**(k - 1)
        return t


class reduced_totient(Function):
    r"""
    Calculate the Carmichael reduced totient function lambda(n)

    ``reduced_totient(n)`` or `\lambda(n)` is the smallest m > 0 such that
    `k^m \equiv 1 \mod n` for all k relatively prime to n.

    Examples
    ========

    >>> from sympy.ntheory import reduced_totient
    >>> reduced_totient(1)
    1
    >>> reduced_totient(8)
    2
    >>> reduced_totient(30)
    4

    See Also
    ========

    totient

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Carmichael_function
    .. [2] https://mathworld.wolfram.com/CarmichaelFunction.html

    """
    @classmethod
    def eval(cls, n):
        if n.is_Integer:
            if n < 1:
                raise ValueError("n must be a positive integer")
            factors = factorint(n)
            return cls._from_factors(factors)

    @classmethod
    def _from_factors(self, factors):
        """Subroutine to compute totient from already-computed factors
        """
        t = 1
        for p, k in factors.items():
            if p == 2 and k > 2:
                t = ilcm(t, 2**(k - 2))
            else:
                t = ilcm(t, (p - 1) * p**(k - 1))
        return t

    @classmethod
    def _from_distinct_primes(self, *args):
        """Subroutine to compute totient from the list of assumed
        distinct primes
        """
        args = [p - 1 for p in args]
        return ilcm(*args)

    def _eval_is_integer(self):
        return fuzzy_and([self.args[0].is_integer, self.args[0].is_positive])


class divisor_sigma(Function):
    r"""
    Calculate the divisor function `\sigma_k(n)` for positive integer n

    ``divisor_sigma(n, k)`` is equal to ``sum([x**k for x in divisors(n)])``

    If n's prime factorization is:

    .. math ::
        n = \prod_{i=1}^\omega p_i^{m_i},

    then

    .. math ::
        \sigma_k(n) = \prod_{i=1}^\omega (1+p_i^k+p_i^{2k}+\cdots
        + p_i^{m_ik}).

    Parameters
    ==========

    n : integer

    k : integer, optional
        power of divisors in the sum

        for k = 0, 1:
        ``divisor_sigma(n, 0)`` is equal to ``divisor_count(n)``
        ``divisor_sigma(n, 1)`` is equal to ``sum(divisors(n))``

        Default for k is 1.

    Examples
    ========

    >>> from sympy.ntheory import divisor_sigma
    >>> divisor_sigma(18, 0)
    6
    >>> divisor_sigma(39, 1)
    56
    >>> divisor_sigma(12, 2)
    210
    >>> divisor_sigma(37)
    38

    See Also
    ========

    divisor_count, totient, divisors, factorint

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Divisor_function

    """

    @classmethod
    def eval(cls, n, k=S.One):
        k = sympify(k)

        if n.is_prime:
            return 1 + n**k

        if n.is_Integer:
            if n <= 0:
                raise ValueError("n must be a positive integer")
            elif k.is_Integer:
                k = int(k)
                return Integer(math.prod(
                    (p**(k*(e + 1)) - 1)//(p**k - 1) if k != 0
                    else e + 1 for p, e in factorint(n).items()))
            else:
                return Mul(*[(p**(k*(e + 1)) - 1)/(p**k - 1) if k != 0
                           else e + 1 for p, e in factorint(n).items()])

        if n.is_integer:  # symbolic case
            args = []
            for p, e in (_.as_base_exp() for _ in Mul.make_args(n)):
                if p.is_prime and e.is_positive:
                    args.append((p**(k*(e + 1)) - 1)/(p**k - 1) if
                                k != 0 else e + 1)
                else:
                    return
            return Mul(*args)


def core(n, t=2):
    r"""
    Calculate core(n, t) = `core_t(n)` of a positive integer n

    ``core_2(n)`` is equal to the squarefree part of n

    If n's prime factorization is:

    .. math ::
        n = \prod_{i=1}^\omega p_i^{m_i},

    then

    .. math ::
        core_t(n) = \prod_{i=1}^\omega p_i^{m_i \mod t}.

    Parameters
    ==========

    n : integer

    t : integer
        core(n, t) calculates the t-th power free part of n

        ``core(n, 2)`` is the squarefree part of ``n``
        ``core(n, 3)`` is the cubefree part of ``n``

        Default for t is 2.

    Examples
    ========

    >>> from sympy.ntheory.factor_ import core
    >>> core(24, 2)
    6
    >>> core(9424, 3)
    1178
    >>> core(379238)
    379238
    >>> core(15**11, 10)
    15

    See Also
    ========

    factorint, sympy.solvers.diophantine.diophantine.square_factor

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Square-free_integer#Squarefree_core

    """

    n = as_int(n)
    t = as_int(t)
    if n <= 0:
        raise ValueError("n must be a positive integer")
    elif t <= 1:
        raise ValueError("t must be >= 2")
    else:
        y = 1
        for p, e in factorint(n).items():
            y *= p**(e % t)
        return y


class udivisor_sigma(Function):
    r"""
    Calculate the unitary divisor function `\sigma_k^*(n)` for positive integer n

    ``udivisor_sigma(n, k)`` is equal to ``sum([x**k for x in udivisors(n)])``

    If n's prime factorization is:

    .. math ::
        n = \prod_{i=1}^\omega p_i^{m_i},

    then

    .. math ::
        \sigma_k^*(n) = \prod_{i=1}^\omega (1+ p_i^{m_ik}).

    Parameters
    ==========

    k : power of divisors in the sum

        for k = 0, 1:
        ``udivisor_sigma(n, 0)`` is equal to ``udivisor_count(n)``
        ``udivisor_sigma(n, 1)`` is equal to ``sum(udivisors(n))``

        Default for k is 1.

    Examples
    ========

    >>> from sympy.ntheory.factor_ import udivisor_sigma
    >>> udivisor_sigma(18, 0)
    4
    >>> udivisor_sigma(74, 1)
    114
    >>> udivisor_sigma(36, 3)
    47450
    >>> udivisor_sigma(111)
    152

    See Also
    ========

    divisor_count, totient, divisors, udivisors, udivisor_count, divisor_sigma,
    factorint

    References
    ==========

    .. [1] https://mathworld.wolfram.com/UnitaryDivisorFunction.html

    """

    @classmethod
    def eval(cls, n, k=S.One):
        k = sympify(k)
        if n.is_prime:
            return 1 + n**k
        if n.is_Integer:
            if n <= 0:
                raise ValueError("n must be a positive integer")
            else:
                return Mul(*[1+p**(k*e) for p, e in factorint(n).items()])


class primenu(Function):
    r"""
    Calculate the number of distinct prime factors for a positive integer n.

    If n's prime factorization is:

    .. math ::
        n = \prod_{i=1}^k p_i^{m_i},

    then ``primenu(n)`` or `\nu(n)` is:

    .. math ::
        \nu(n) = k.

    Examples
    ========

    >>> from sympy.ntheory.factor_ import primenu
    >>> primenu(1)
    0
    >>> primenu(30)
    3

    See Also
    ========

    factorint

    References
    ==========

    .. [1] https://mathworld.wolfram.com/PrimeFactor.html

    """

    @classmethod
    def eval(cls, n):
        if n.is_Integer:
            if n <= 0:
                raise ValueError("n must be a positive integer")
            else:
                return len(factorint(n).keys())


class primeomega(Function):
    r"""
    Calculate the number of prime factors counting multiplicities for a
    positive integer n.

    If n's prime factorization is:

    .. math ::
        n = \prod_{i=1}^k p_i^{m_i},

    then ``primeomega(n)``  or `\Omega(n)` is:

    .. math ::
        \Omega(n) = \sum_{i=1}^k m_i.

    Examples
    ========

    >>> from sympy.ntheory.factor_ import primeomega
    >>> primeomega(1)
    0
    >>> primeomega(20)
    3

    See Also
    ========

    factorint

    References
    ==========

    .. [1] https://mathworld.wolfram.com/PrimeFactor.html

    """

    @classmethod
    def eval(cls, n):
        if n.is_Integer:
            if n <= 0:
                raise ValueError("n must be a positive integer")
            else:
                return sum(factorint(n).values())


def mersenne_prime_exponent(nth):
    """Returns the exponent ``i`` for the nth Mersenne prime (which
    has the form `2^i - 1`).

    Examples
    ========

    >>> from sympy.ntheory.factor_ import mersenne_prime_exponent
    >>> mersenne_prime_exponent(1)
    2
    >>> mersenne_prime_exponent(20)
    4423
    """
    n = as_int(nth)
    if n < 1:
        raise ValueError("nth must be a positive integer; mersenne_prime_exponent(1) == 2")
    if n > 51:
        raise ValueError("There are only 51 perfect numbers; nth must be less than or equal to 51")
    return MERSENNE_PRIME_EXPONENTS[n - 1]


def is_perfect(n):
    """Returns True if ``n`` is a perfect number, else False.

    A perfect number is equal to the sum of its positive, proper divisors.

    Examples
    ========

    >>> from sympy.ntheory.factor_ import is_perfect, divisors, divisor_sigma
    >>> is_perfect(20)
    False
    >>> is_perfect(6)
    True
    >>> 6 == divisor_sigma(6) - 6 == sum(divisors(6)[:-1])
    True

    References
    ==========

    .. [1] https://mathworld.wolfram.com/PerfectNumber.html
    .. [2] https://en.wikipedia.org/wiki/Perfect_number

    """

    n = as_int(n)
    if _isperfect(n):
        return True

    # all perfect numbers for Mersenne primes with exponents
    # less than or equal to 43112609 are known
    iknow = MERSENNE_PRIME_EXPONENTS.index(43112609)
    if iknow <= len(PERFECT) - 1 and n <= PERFECT[iknow]:
        # there may be gaps between this and larger known values
        # so only conclude in the range for which all values
        # are known
        return False
    if n%2 == 0:
        last2 = n % 100
        if last2 != 28 and last2 % 10 != 6:
            return False
        r, b = integer_nthroot(1 + 8*n, 2)
        if not b:
            return False
        m, x = divmod(1 + r, 4)
        if x:
            return False
        e, b = integer_log(m, 2)
        if not b:
            return False
    else:
        if n < 10**2000:  # https://www.lirmm.fr/~ochem/opn/
            return False
        if n % 105 == 0:  # not divis by 105
            return False
        if not any(n%m == r for m, r in [(12, 1), (468, 117), (324, 81)]):
            return False
        # there are many criteria that the factor structure of n
        # must meet; since we will have to factor it to test the
        # structure we will have the factors and can then check
        # to see whether it is a perfect number or not. So we
        # skip the structure checks and go straight to the final
        # test below.
    rv = divisor_sigma(n) - n
    if rv == n:
        if n%2 == 0:
            raise ValueError(filldedent('''
                This even number is perfect and is associated with a
                Mersenne Prime, 2^%s - 1. It should be
                added to SymPy.''' % (e + 1)))
        else:
            raise ValueError(filldedent('''In 1888, Sylvester stated: "
                ...a prolonged meditation on the subject has satisfied
                me that the existence of any one such [odd perfect number]
                -- its escape, so to say, from the complex web of conditions
                which hem it in on all sides -- would be little short of a
                miracle." I guess SymPy just found that miracle and it
                factors like this: %s''' % factorint(n)))


def is_mersenne_prime(n):
    """Returns True if  ``n`` is a Mersenne prime, else False.

    A Mersenne prime is a prime number having the form `2^i - 1`.

    Examples
    ========

    >>> from sympy.ntheory.factor_ import is_mersenne_prime
    >>> is_mersenne_prime(6)
    False
    >>> is_mersenne_prime(127)
    True

    References
    ==========

    .. [1] https://mathworld.wolfram.com/MersennePrime.html

    """

    n = as_int(n)
    if _ismersenneprime(n):
        return True
    if not isprime(n):
        return False
    r, b = integer_log(n + 1, 2)
    if not b:
        return False
    raise ValueError(filldedent('''
        This Mersenne Prime, 2^%s - 1, should
        be added to SymPy's known values.''' % r))


def abundance(n):
    """Returns the difference between the sum of the positive
    proper divisors of a number and the number.

    Examples
    ========

    >>> from sympy.ntheory import abundance, is_perfect, is_abundant
    >>> abundance(6)
    0
    >>> is_perfect(6)
    True
    >>> abundance(10)
    -2
    >>> is_abundant(10)
    False
    """
    return divisor_sigma(n, 1) - 2 * n


def is_abundant(n):
    """Returns True if ``n`` is an abundant number, else False.

    A abundant number is smaller than the sum of its positive proper divisors.

    Examples
    ========

    >>> from sympy.ntheory.factor_ import is_abundant
    >>> is_abundant(20)
    True
    >>> is_abundant(15)
    False

    References
    ==========

    .. [1] https://mathworld.wolfram.com/AbundantNumber.html

    """
    n = as_int(n)
    if is_perfect(n):
        return False
    return n % 6 == 0 or bool(abundance(n) > 0)


def is_deficient(n):
    """Returns True if ``n`` is a deficient number, else False.

    A deficient number is greater than the sum of its positive proper divisors.

    Examples
    ========

    >>> from sympy.ntheory.factor_ import is_deficient
    >>> is_deficient(20)
    False
    >>> is_deficient(15)
    True

    References
    ==========

    .. [1] https://mathworld.wolfram.com/DeficientNumber.html

    """
    n = as_int(n)
    if is_perfect(n):
        return False
    return bool(abundance(n) < 0)


def is_amicable(m, n):
    """Returns True if the numbers `m` and `n` are "amicable", else False.

    Amicable numbers are two different numbers so related that the sum
    of the proper divisors of each is equal to that of the other.

    Examples
    ========

    >>> from sympy.ntheory.factor_ import is_amicable, divisor_sigma
    >>> is_amicable(220, 284)
    True
    >>> divisor_sigma(220) == divisor_sigma(284)
    True

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Amicable_numbers

    """
    if m == n:
        return False
    a, b = (divisor_sigma(i) for i in (m, n))
    return a == b == (m + n)


def dra(n, b):
    """
    Returns the additive digital root of a natural number ``n`` in base ``b``
    which is a single digit value obtained by an iterative process of summing
    digits, on each iteration using the result from the previous iteration to
    compute a digit sum.

    Examples
    ========

    >>> from sympy.ntheory.factor_ import dra
    >>> dra(3110, 12)
    8

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Digital_root

    """

    num = abs(as_int(n))
    b = as_int(b)
    if b <= 1:
        raise ValueError("Base should be an integer greater than 1")

    if num == 0:
        return 0

    return (1 + (num - 1) % (b - 1))


def drm(n, b):
    """
    Returns the multiplicative digital root of a natural number ``n`` in a given
    base ``b`` which is a single digit value obtained by an iterative process of
    multiplying digits, on each iteration using the result from the previous
    iteration to compute the digit multiplication.

    Examples
    ========

    >>> from sympy.ntheory.factor_ import drm
    >>> drm(9876, 10)
    0

    >>> drm(49, 10)
    8

    References
    ==========

    .. [1] https://mathworld.wolfram.com/MultiplicativeDigitalRoot.html

    """

    n = abs(as_int(n))
    b = as_int(b)
    if b <= 1:
        raise ValueError("Base should be an integer greater than 1")
    while n > b:
        mul = 1
        while n > 1:
            n, r = divmod(n, b)
            if r == 0:
                return 0
            mul *= r
        n = mul
    return n
