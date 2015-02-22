"""
Primality testing

"""

from __future__ import print_function, division
from sympy.core.compatibility import range

# pseudoprimes that will pass through last mr_safe test
_pseudos = set([
            669094855201,
           1052516956501, 2007193456621, 2744715551581, 9542968210729,
          17699592963781, 19671510288601,
          24983920772821, 24984938689453, 29661584268781, 37473222618541,
          46856248255981, 47922612926653, 48103703944453, 49110566041153,
          49752242681221, 91206655032481, 91481980096033, 119034193492321,
         123645258399601, 128928036060253, 137364148720147, 150753857310253,
         153131886327421, 155216912613121, 185610214763821, 224334357392701,
         227752294950181, 230058334559041, 304562854940401, 306001576998253,
         335788261073821, 377133492079081, 379242177424951, 389970770948461,
         397319638319521, 448114903362253, 523235160050221, 628999496281621,
         699349238838253, 746667678235753, 790198268451301, 794036495175661,
         823820871230281, 867739535711821, 1039918661294761, 1099127938585141,
        1104388025338153, 1173374598605653, 1262797719066157, 1265872947674653,
        1325898212229667, 1327034517143653, 1418575746675583, 1666122072463621,
        1837400535259453, 1857422490084961, 1870756820971741, 1914550540480717,
        2018963273468221, 2163829000939453, 2206020317369221, 2301037384029121,
        2416062055125421, 2435076500074921, 2545656135020833, 2594428516569781,
        2669983768115821, 2690937050990653, 2758640869506607, 2833525461416653,
        2876662942007221, 2932155806957821, 2957010595723801, 3183606449929153,
        3220133449185901, 3424103775720253, 3625360152399541, 3939300299037421,
        3947917710714841, 3980273496750253, 4182256679324041, 4450605887818261,
        4727893739521501, 4750350311306953, 4755334362931153, 5756440863559753,
        5760976603475341, 5794399356078761, 5954850603819253, 6125544931991761,
        6320931714094861, 6347593619672581, 6406268028524101, 6510632945054941,
        6620082224794741, 6627325072566061, 6844056606431101, 6989404981060153,
        7144293947609521, 7288348593229021, 7288539837129253, 7406102904971689,
        7430233301822341, 7576425305871193, 7601696719033861, 7803926845356487,
        7892007967006633, 7947797946559453, 8207000460596953, 8295064717807513,
        8337196000698841, 8352714234009421, 8389755717406381, 8509654470665701,
        8757647355282841, 8903933671696381, 8996133652295653, 9074421465661261,
        9157536631454221, 9188353522314541])


def _test(n, base, s, t):
    """Miller-Rabin strong pseudoprime test for one base.
    Return False if n is definitely composite, True if n is
    probably prime, with a probability greater than 3/4.

    """
    # do the Fermat test
    b = pow(base, t, n)
    if b == 1 or b == n - 1:
        return True
    else:
        for j in range(1, s):
            b = pow(b, 2, n)
            if b == n - 1:
                return True
            # see I. Niven et al. "An Introduction to Theory of Numbers", page 78
            if b == 1:
                return False
    return False


def mr(n, bases):
    """Perform a Miller-Rabin strong pseudoprime test on n using a
    given list of bases/witnesses.

    References
    ==========

    - Richard Crandall & Carl Pomerance (2005), "Prime Numbers:
      A Computational Perspective", Springer, 2nd edition, 135-138

    A list of thresholds and the bases they require are here:
    http://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Deterministic_variants_of_the_test

    Examples
    ========

    >>> from sympy.ntheory.primetest import mr
    >>> mr(1373651, [2, 3])
    False
    >>> mr(479001599, [31, 73])
    True

    """
    from sympy.ntheory.factor_ import trailing
    from sympy.polys.domains import ZZ

    n = int(n)
    if n < 2:
        return False
    # remove powers of 2 from n = t * 2**s + 1
    s = trailing(n - 1)
    t = n >> s
    for base in bases:
        base = ZZ(base)
        if not _test(n, base, s, t):
            return False
    return True


def _mr_safe(n):
    """For n < 10**16, use the Miller-Rabin test to determine with
    certainty (unless the code is buggy!) whether n is prime.

    Although the primes 2 through 17 are sufficient to confirm that a number
    less than 341550071728322 (that is not prime 2 through 17) is prime, this
    range is broken up into smaller ranges with earlier ranges requiring less
    work. For example, for n < 1373653 only the bases 2 and 3 need be tested.

    What makes this a "safe" Miller-Rabin routine is that for n less than
    the indicated limit, the given bases have been confirmed to detect all
    composite numbers. What can potentially make this routine "unsafe" is
    including ranges for which previous tests do not removes prime factors of
    the bases being used. For example, this routine assumes that 2 and 3 have
    already been removed as prime; but if the first test were the one for
    n < 170584961 (that uses bases 350 and 3958281543) the routine would have
    to ensure that the primes 5, 7, 29, 67, 679067 are already removed or else
    they will be reported as being composite. For this reason it is helpful to
    list the prime factors of the bases being tested as is done below. The
    _mr_safe_helper can be used to generate this info-tag.

    References for the bounds:
    ==========================

    1. http://primes.utm.edu/prove/prove2_3.html
    2. http://www.trnicely.net/misc/mpzspsp.html
    3. http://en.wikipedia.org/wiki/Miller-Rabin_primality_test#
        Accuracy_of_the_test
    4. http://primes.utm.edu/glossary/xpage/Pseudoprime.html
    5. http://uucode.com/obf/dalbec/alg.html#sprp

    """

    if n < 1373653:
        return mr(n, [2, 3])
        #[2, 3] stot = 1 clear == bases
        # these two (and similar below) are commented out since they are
        # more expensive in terms of stot than a later test.
        #if n < 9080191: return mr(n, [31, 73]) # ref [3]
        # [31, 73] stot = 4 clear == bases
        #if n < 25326001: return mr(n, [2, 3, 5])
        # [2, 3, 5] stot = 3 clear == bases
    if n < 170584961:
        return mr(n, [350, 3958281543])
        # [350, 3958281543] stot = 1 clear [2, 3, 5, 7, 29, 67, 679067]
    if n < 4759123141:
        return mr(n, [2, 7, 61])  # ref [3]
        # [2, 7, 61] stot = 3 clear == bases
    if n < 75792980677:
        return mr(n, [2, 379215, 457083754])
        # [2, 379215, 457083754] stot = 1 clear [2, 3, 5, 53, 228541877]
        #if n < 118670087467: return n is not 3215031751 and mr(n, [2, 3, 5, 7]) # ref [3]
        # [2, 3, 5, 7] stot = 4 clear == bases
    if n < 1000000000000:
        return mr(n, [2, 13, 23, 1662803])
        # [2, 13, 23, 1662803] stot = 4 clear == bases
        #if n < 2152302898747: return mr(n, [2, 3, 5, 7, 11])
        # [2, 3, 5, 7, 11] stot = 5 clear == bases
        #if n < 3474749660383: return mr(n, [2, 3, 5, 7, 11, 13])
        # [2, 3, 5, 7, 11, 13] stot = 7 clear == bases
        #if n < 21652684502221: return mr(n, [2, 1215, 34862, 574237825])
        # [2, 1215, 34862, 574237825] stot = 8 clear [2, 3, 5, 7, 17431, 3281359]
        #if n < 341550071728321: return mr(n, [2, 3, 5, 7, 11, 13, 17])
        # [2, 3, 5, 7, 11, 13, 17] stot = 11 clear == bases
    if n < 10000000000000000:
        return mr(n, [2, 3, 7, 61, 24251]) and n not in _pseudos
        # [2, 3, 7, 61, 24251] stot = 5 clear == bases
    raise ValueError("n too large")


def isprime(n):
    """
    Test if n is a prime number (True) or not (False). For n < 10**16 the
    answer is accurate; greater n values have a small probability of actually
    being pseudoprimes.

    Negative primes (e.g. -2) are not considered prime.

    The function first looks for trivial factors, and if none is found,
    performs a safe Miller-Rabin strong pseudoprime test with bases
    that are known to prove a number prime. Finally, a general Miller-Rabin
    test is done with the first k bases which will report a pseudoprime as a
    prime with an error of about 4**-k. The current value of k is 46 so the
    error is about 2 x 10**-28.

    Examples
    ========

    >>> from sympy.ntheory import isprime
    >>> isprime(13)
    True
    >>> isprime(15)
    False

    See Also
    ========

    sympy.ntheory.generate.primerange : Generates all primes in a given range
    sympy.ntheory.generate.primepi : Return the number of primes less than or equal to n
    sympy.ntheory.generate.prime : Return the nth prime

    """
    n = int(n)
    if n < 2:
        return False
    if n & 1 == 0:
        return n == 2
    if n <= 23001:
        return pow(2, n, n) == 2 and n not in [341, 561, 645, 1105, 1387, 1729,
                                               1905, 2047, 2465, 2701, 2821,
                                               3277, 4033, 4369, 4371, 4681,
                                               5461, 6601, 7957, 8321, 8481,
                                               8911, 10261, 10585, 11305,
                                               12801, 13741, 13747, 13981,
                                               14491, 15709, 15841, 16705,
                                               18705, 18721, 19951, 23001]
    try:
        return _mr_safe(n)
    except ValueError:
        # prime list to use when number must be tested as a probable prime;
        # these are the 46 primes less than 200
        bases = [
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
            53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
            109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
            173, 179, 181, 191, 193, 197, 199]
        return mr(n, bases)


def _mr_safe_helper(_s):
    """
    Analyze a (new) mr_safe line for for total number of s's to
    be tested in _test along with the primes that must be cleared
    by a previous test.

    e.g.
    >>> from sympy.ntheory.primetest import _mr_safe_helper
    >>> print(_mr_safe_helper("if n < 170584961: return mr(n, [350, 3958281543])"))
     # [350, 3958281543] stot = 1 clear [2, 3, 5, 7, 29, 67, 679067]
    >>> print(_mr_safe_helper('return mr(n, [2, 379215, 457083754])'))
     # [2, 379215, 457083754] stot = 1 clear [2, 3, 5, 53, 228541877]

    """

    def _info(bases):
        """
        Analyze the list of bases, reporting the number of 'j-loops' that
        will be required if this list is passed to _test (stot) and the primes
        that must be cleared by a previous test.

        This info tag should then be appended to any new mr_safe line
        that is added so someone can easily see whether that line satisfies
        the requirements of mr_safe (see docstring there for details).

        """
        from sympy.ntheory.factor_ import factorint, trailing

        factors = []
        tot = 0
        for b in bases:
            tot += trailing(b - 1)
            f = factorint(b)
            factors.extend(f)
        factors = sorted(set(factors))
        bases = sorted(set(bases))
        if bases == factors:
            factors = '== bases'
        else:
            factors = str(factors)
        return ' # %s stot = %s clear %s' % tuple(
            [str(x).replace('L', '') for x in (list(bases), tot, factors)])

    _r = [int(_x) for _x in _s.split('[')[1].split(']')[0].split(',')]
    return _info(_r)
