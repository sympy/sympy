"""
Primality testing

"""

_tiny_primes = [2, 3, 5, 7, 11, 13, 17, 19]
_max_tiny_prime = 19
_tiny_primes_set = set(_tiny_primes)

# prime list to use when number must be tested as a probable prime.
#>>> list(primerange(2, 200))
_isprime_fallback_primes = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
    53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
    109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
    173, 179, 181, 191, 193, 197, 199]
#>>> len(_)
#46
# pseudoprimes that will pass through last mr_safe test
pseudos = set([
        669094855201,
        1052516956501,2007193456621,2744715551581,9542968210729,
        17699592963781,19671510288601,
        24983920772821,24984938689453,29661584268781,37473222618541,
        46856248255981,47922612926653,48103703944453,49110566041153,
        49752242681221,91206655032481,91481980096033,119034193492321,
        123645258399601,128928036060253,137364148720147,150753857310253,
        153131886327421,155216912613121,185610214763821,224334357392701,
        227752294950181,230058334559041,304562854940401,306001576998253,
        335788261073821,377133492079081,379242177424951,389970770948461,
        397319638319521,448114903362253,523235160050221,628999496281621,
        699349238838253,746667678235753,790198268451301,794036495175661,
        823820871230281,867739535711821,1039918661294761,1099127938585141,
        1104388025338153,1173374598605653,1262797719066157,1265872947674653,
        1325898212229667,1327034517143653,1418575746675583,1666122072463621,
        1837400535259453,1857422490084961,1870756820971741,1914550540480717,
        2018963273468221,2163829000939453,2206020317369221,2301037384029121,
        2416062055125421,2435076500074921,2545656135020833,2594428516569781,
        2669983768115821,2690937050990653,2758640869506607,2833525461416653,
        2876662942007221,2932155806957821,2957010595723801,3183606449929153,
        3220133449185901,3424103775720253,3625360152399541,3939300299037421,
        3947917710714841,3980273496750253,4182256679324041,4450605887818261,
        4727893739521501,4750350311306953,4755334362931153,5756440863559753,
        5760976603475341,5794399356078761,5954850603819253,6125544931991761,
        6320931714094861,6347593619672581,6406268028524101,6510632945054941,
        6620082224794741,6627325072566061,6844056606431101,6989404981060153,
        7144293947609521,7288348593229021,7288539837129253,7406102904971689,
        7430233301822341,7576425305871193,7601696719033861,7803926845356487,
        7892007967006633,7947797946559453,8207000460596953,8295064717807513,
        8337196000698841,8352714234009421,8389755717406381,8509654470665701,
        8757647355282841,8903933671696381,8996133652295653,9074421465661261,
        9157536631454221,9188353522314541])

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
    probably prime, with a probability greater than 3/4.
    """

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

def _mr_safe(n):
    """For n < 1e16, use the Miller-Rabin test to determine with
    certainty (unless the code is buggy!) whether n is prime.

    NOTE: If the list of bases to test includes a base that is composite,
    then n must not contain any factors of the bases or
    else an incorrect result will be returned. e.g. 1943 (the product of
    67 and 29) and the prime 679067 cannot be tested with the bases 350 and
    3958281543 (specifically, the latter) since the factors of the latter are
    29, 67 and 679067. The way to make sure that a valid n is tested is
    to make sure that a preceding test allows only n values greater than
    the largest prime factor in a successive test e.g. testing n < 1373653
    before testing with bases 350 and 3958281543 guarantees that n will not
    have any factor less than 1373653.

    If new ranges are added below it would be nice to indicate what the largest
    prime of a given base set is so it is clear whether this condition is being
    met.

    Reference for the bounds:
    [1] http://primes.utm.edu/prove/prove2_3.html
    [2] http://www.trnicely.net/misc/mpzspsp.html
    [3] http://en.wikipedia.org/wiki/Miller-Rabin_primality_test#
        Accuracy_of_the_test
    [4] http://zdu.spaces.live.com/?_c11_BlogPart_pagedir=
        Next&_c11_BlogPart_handle=cns!C95152CB25EF2037!
        138&_c11_BlogPart_BlogPart=blogview&_c=BlogPart
    [5] http://primes.utm.edu/glossary/xpage/Pseudoprime.html
    [6] http://uucode.com/obf/dalbec/alg.html#sprp
    """

    n = int(n)
    if n < 1373653: return mr(n, [2, 3])
    if n < 170584961: return mr(n, [350, 3958281543])
    if n < 4759123141: return mr(n, [2, 7, 61]) # ref [3]
    if n < 75792980677: return mr(n, [2, 379215, 457083754])
    if n < 1000000000000: return mr(n, [2, 13, 23, 1662803])
    if n < 10000000000000000: return mr(n, [2, 3, 7, 61, 24251]) \
       and n not in pseudos
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
    test is done with the first k bases which, which will report a
    pseudoprime as a prime with an error of about 4**-k. The current value
    of k is 46 so the error is about 2 x 10**-28.

    Example usage
    =============
    >>> from sympy.ntheory import isprime
    >>> isprime(13)
        True
        >>> isprime(15)
        False

    """
    n = int(n)
    if n < 2:
        return False
    if n & 1 == 0:
        return n == 2
    if _is_tiny_prime(n):
        return True
    if _has_tiny_factor(n):
        return False
    try:
        return _mr_safe(n)
    except ValueError:
        return mr(n, _isprime_fallback_primes)

def _mr_safe_helper(_s):
    """
    Analyze a (new) mr_safe line for for total number of s's to
    be tested in _test along with the highest prime that occurs
    in the list within that line.

    e.g.
    >>> from sympy.ntheory.primetest import _mr_safe_helper
    >>> print _mr_safe_helper("if n < 170584961: return mr(n, [350, 3958281543])")
     # [350, 3958281543] stot = 1 pmax = 1319427181
    """

    def _maxfac(smalln):
        """
        Find the largest prime factor of smalln without resorting to
        sympy's own functions.
        """
        from math import sqrt
        assert smalln > 1
        assert _max_tiny_prime > 2
        maxp = int(sqrt(smalln))
        for p in _tiny_primes:
            if p > maxp:
                return smalln
            if smalln % p == 0:
                return max([smalln//p, p])
        for p in range(_max_tiny_prime + 2, maxp, 2):
            if smalln % p == 0:
                return max([smalln//p, p])
        return smalln

    def _info(bases):
        """
        Analyze the list of bases, reporting the number of 'j-loops' that
        will be required if this list is passed to _test and the maximum
        prime that occurs (perhaps as a factor of a composite base) in the
        list. This info tag should then be appended to any new mr_safe line
        that is added so someone can easily see whether that line satisfies
        the requirements of mr_safe (see docstring there for details).
        """
        pmax = 0
        tot = 0
        for b in bases:
            tot += _factor_pow2(b-1)[0]
            pmax = max([pmax, _maxfac(b)])
        return ' # %s stot = %s pmax = %s' % tuple(
                [str(x).replace('L','') for x in (list(bases), tot, pmax)])

    _r = [int(_x) for _x in _s.split('[')[1].split(']')[0].split(',')]
    return _info(_r)
