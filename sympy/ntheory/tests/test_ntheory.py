from collections import defaultdict
from sympy import Sieve, binomial_coefficients, binomial_coefficients_list, \
    multinomial_coefficients, Mul, S, Pow, sieve, Symbol, summation, \
    factorial as fac, pi, GoldenRatio as phi, sqrt
from sympy.core.numbers import Integer, Rational
from sympy.core.compatibility import long, range

from sympy.ntheory import isprime, n_order, is_primitive_root, \
    is_quad_residue, legendre_symbol, jacobi_symbol, npartitions, totient, \
    factorint, primefactors, divisors, randprime, nextprime, prevprime, \
    primerange, primepi, prime, pollard_rho, perfect_power, multiplicity, \
    trailing, divisor_count, primorial, pollard_pm1, \
    sqrt_mod, primitive_root, quadratic_residues, is_nthpow_residue, \
    nthroot_mod, sqrt_mod_iter, mobius, divisor_sigma

from sympy.ntheory.residue_ntheory import _primitive_root_prime_iter
from sympy.ntheory.factor_ import smoothness, smoothness_p, \
    antidivisors, antidivisor_count, core
from sympy.ntheory.generate import cycle_length
from sympy.ntheory.primetest import _mr_safe_helper, mr
from sympy.ntheory.bbp_pi import pi_hex_digits
from sympy.ntheory.modular import crt, crt1, crt2, solve_congruence
from sympy.ntheory.continued_fraction import \
    (continued_fraction_periodic as cf_p,
     continued_fraction_iterator as cf_i,
     continued_fraction_convergents as cf_c,
     continued_fraction_reduce as cf_r)
from sympy.ntheory.egyptian_fraction import egyptian_fraction

from sympy.core.add import Add

from sympy.polys.domains import ZZ

from sympy.utilities.pytest import raises
from sympy.utilities.iterables import capture
from sympy.utilities.randtest import random_complex_number
from sympy.ntheory.multinomial import multinomial_coefficients_iterator

def test_trailing():
    assert trailing(0) == 0
    assert trailing(1) == 0
    assert trailing(-1) == 0
    assert trailing(2) == 1
    assert trailing(7) == 0
    assert trailing(-7) == 0
    for i in range(100):
        assert trailing((1 << i)) == i
        assert trailing((1 << i) * 31337) == i
    assert trailing((1 << 1000001)) == 1000001
    assert trailing((1 << 273956)*7**37) == 273956


def test_multiplicity():
    for b in range(2, 20):
        for i in range(100):
            assert multiplicity(b, b**i) == i
            assert multiplicity(b, (b**i) * 23) == i
            assert multiplicity(b, (b**i) * 1000249) == i
    # Should be fast
    assert multiplicity(10, 10**10023) == 10023
    # Should exit quickly
    assert multiplicity(10**10, 10**10) == 1
    # Should raise errors for bad input
    raises(ValueError, lambda: multiplicity(1, 1))
    raises(ValueError, lambda: multiplicity(1, 2))
    raises(ValueError, lambda: multiplicity(1.3, 2))
    raises(ValueError, lambda: multiplicity(5, 0))

    # handles Rationals
    assert multiplicity(10, Rational(30, 7)) == 0
    assert multiplicity(Rational(2, 7), Rational(4, 7)) == 1
    assert multiplicity(Rational(1, 7), Rational(3, 49)) == 2
    assert multiplicity(Rational(2, 7), Rational(7, 2)) == -1
    assert multiplicity(3, Rational(1, 9)) == -2


def test_perfect_power():
    assert perfect_power(0) is False
    assert perfect_power(1) is False
    assert perfect_power(2) is False
    assert perfect_power(3) is False
    assert perfect_power(4) == (2, 2)
    assert perfect_power(14) is False
    assert perfect_power(25) == (5, 2)
    assert perfect_power(22) is False
    assert perfect_power(22, [2]) is False
    assert perfect_power(137**(3*5*13)) == (137, 3*5*13)
    assert perfect_power(137**(3*5*13) + 1) is False
    assert perfect_power(137**(3*5*13) - 1) is False
    assert perfect_power(103005006004**7) == (103005006004, 7)
    assert perfect_power(103005006004**7 + 1) is False
    assert perfect_power(103005006004**7 - 1) is False
    assert perfect_power(103005006004**12) == (103005006004, 12)
    assert perfect_power(103005006004**12 + 1) is False
    assert perfect_power(103005006004**12 - 1) is False
    assert perfect_power(2**10007) == (2, 10007)
    assert perfect_power(2**10007 + 1) is False
    assert perfect_power(2**10007 - 1) is False
    assert perfect_power((9**99 + 1)**60) == (9**99 + 1, 60)
    assert perfect_power((9**99 + 1)**60 + 1) is False
    assert perfect_power((9**99 + 1)**60 - 1) is False
    assert perfect_power((10**40000)**2, big=False) == (10**40000, 2)
    assert perfect_power(10**100000) == (10, 100000)
    assert perfect_power(10**100001) == (10, 100001)
    assert perfect_power(13**4, [3, 5]) is False
    assert perfect_power(3**4, [3, 10], factor=0) is False
    assert perfect_power(3**3*5**3) == (15, 3)
    assert perfect_power(2**3*5**5) is False
    assert perfect_power(2*13**4) is False
    assert perfect_power(2**5*3**3) is False


def test_isprime():
    s = Sieve()
    s.extend(100000)
    ps = set(s.primerange(2, 100001))
    for n in range(100001):
        # if (n in ps) != isprime(n): print n
        assert (n in ps) == isprime(n)
    assert isprime(179424673)
    # Some Mersenne primes
    assert isprime(2**61 - 1)
    assert isprime(2**89 - 1)
    assert isprime(2**607 - 1)
    assert not isprime(2**601 - 1)
    #Arnault's number
    assert isprime(int('''
803837457453639491257079614341942108138837688287558145837488917522297\
427376533365218650233616396004545791504202360320876656996676098728404\
396540823292873879185086916685732826776177102938969773947016708230428\
687109997439976544144845341155872450633409279022275296229414984230688\
1685404326457534018329786111298960644845216191652872597534901'''))
    # pseudoprime that passes the base set [2, 3, 7, 61, 24251]
    assert not isprime(9188353522314541)

    assert _mr_safe_helper(
        "if n < 170584961: return mr(n, [350, 3958281543])") == \
        ' # [350, 3958281543] stot = 1 clear [2, 3, 5, 7, 29, 67, 679067]'
    assert _mr_safe_helper(
        "if n < 3474749660383: return mr(n, [2, 3, 5, 7, 11, 13])") == \
        ' # [2, 3, 5, 7, 11, 13] stot = 7 clear == bases'


def test_prime():
    assert prime(1) == 2
    assert prime(2) == 3
    assert prime(5) == 11
    assert prime(11) == 31
    assert prime(57) == 269
    assert prime(296) == 1949
    assert prime(559) == 4051
    assert prime(3000) == 27449
    assert prime(4096) == 38873
    assert prime(9096) == 94321
    assert prime(25023) == 287341
    raises(ValueError, lambda: prime(0))


def test_primepi():
    assert primepi(1) == 0
    assert primepi(2) == 1
    assert primepi(5) == 3
    assert primepi(11) == 5
    assert primepi(57) == 16
    assert primepi(296) == 62
    assert primepi(559) == 102
    assert primepi(3000) == 430
    assert primepi(4096) == 564
    assert primepi(9096) == 1128
    assert primepi(25023) == 2763


def test_generate():
    assert nextprime(-4) == 2
    assert nextprime(2) == 3
    assert nextprime(5) == 7
    assert nextprime(12) == 13
    assert nextprime(90) == 97
    assert nextprime(10**40) == (10**40 + 121)
    assert prevprime(3) == 2
    assert prevprime(7) == 5
    assert prevprime(13) == 11
    assert prevprime(97) == 89
    assert prevprime(10**40) == (10**40 - 17)
    assert list(primerange(2, 7)) == [2, 3, 5]
    assert list(primerange(2, 10)) == [2, 3, 5, 7]
    assert list(primerange(1050, 1100)) == [1051, 1061,
        1063, 1069, 1087, 1091, 1093, 1097]
    s = Sieve()
    for i in range(30, 2350, 376):
        for j in range(2, 5096, 1139):
            A = list(s.primerange(i, i + j))
            B = list(primerange(i, i + j))
            assert A == B
    s = Sieve()
    assert s[10] == 29

    assert nextprime(2, 2) == 5

    raises(ValueError, lambda: totient(0))

    raises(ValueError, lambda: primorial(0))

    assert mr(1, [2]) is False

    func = lambda i: (i**2 + 1) % 51
    assert next(cycle_length(func, 4)) == (6, 2)
    assert list(cycle_length(func, 4, values=True)) == \
        [17, 35, 2, 5, 26, 14, 44, 50, 2, 5, 26, 14]
    assert next(cycle_length(func, 4, nmax=5)) == (5, None)
    assert list(cycle_length(func, 4, nmax=5, values=True)) == \
        [17, 35, 2, 5, 26]


def test_randprime():
    import random
    random.seed(1234)
    assert randprime(2, 3) == 2
    assert randprime(1, 3) == 2
    assert randprime(3, 5) == 3
    raises(ValueError, lambda: randprime(20, 22))
    for a in [100, 300, 500, 250000]:
        for b in [100, 300, 500, 250000]:
            p = randprime(a, a + b)
            assert a <= p < (a + b) and isprime(p)


def fac_multiplicity(n, p):
    """Return the power of the prime number p in the
    factorization of n!"""
    if p > n:
        return 0
    if p > n//2:
        return 1
    q, m = n, 0
    while q >= p:
        q //= p
        m += q
    return m


def multiproduct(seq=(), start=1):
    """
    Return the product of a sequence of factors with multiplicities,
    times the value of the parameter ``start``. The input may be a
    sequence of (factor, exponent) pairs or a dict of such pairs.

        >>> multiproduct({3:7, 2:5}, 4) # = 3**7 * 2**5 * 4
        279936

    """
    if not seq:
        return start
    if isinstance(seq, dict):
        seq = iter(seq.items())
    units = start
    multi = []
    for base, exp in seq:
        if not exp:
            continue
        elif exp == 1:
            units *= base
        else:
            if exp % 2:
                units *= base
            multi.append((base, exp//2))
    return units * multiproduct(multi)**2


def test_factorint():
    assert primefactors(123456) == [2, 3, 643]
    assert factorint(0) == {0: 1}
    assert factorint(1) == {}
    assert factorint(-1) == {-1: 1}
    assert factorint(-2) == {-1: 1, 2: 1}
    assert factorint(-16) == {-1: 1, 2: 4}
    assert factorint(2) == {2: 1}
    assert factorint(126) == {2: 1, 3: 2, 7: 1}
    assert factorint(123456) == {2: 6, 3: 1, 643: 1}
    assert factorint(5951757) == {3: 1, 7: 1, 29: 2, 337: 1}
    assert factorint(64015937) == {7993: 1, 8009: 1}
    assert factorint(2**(2**6) + 1) == {274177: 1, 67280421310721: 1}
    assert multiproduct(factorint(fac(200))) == fac(200)
    for b, e in factorint(fac(150)).items():
        assert e == fac_multiplicity(150, b)
    assert factorint(103005006059**7) == {103005006059: 7}
    assert factorint(31337**191) == {31337: 191}
    assert factorint(2**1000 * 3**500 * 257**127 * 383**60) == \
        {2: 1000, 3: 500, 257: 127, 383: 60}
    assert len(factorint(fac(10000))) == 1229
    assert factorint(12932983746293756928584532764589230) == \
        {2: 1, 5: 1, 73: 1, 727719592270351: 1, 63564265087747: 1, 383: 1}
    assert factorint(727719592270351) == {727719592270351: 1}
    assert factorint(2**64 + 1, use_trial=False) == factorint(2**64 + 1)
    for n in range(60000):
        assert multiproduct(factorint(n)) == n
    assert pollard_rho(2**64 + 1, seed=1) == 274177
    assert pollard_rho(19, seed=1) is None
    assert factorint(3, limit=2) == {3: 1}
    assert factorint(12345) == {3: 1, 5: 1, 823: 1}
    assert factorint(
        12345, limit=3) == {4115: 1, 3: 1}  # the 5 is greater than the limit
    assert factorint(1, limit=1) == {}
    assert factorint(0, 3) == {0: 1}
    assert factorint(12, limit=1) == {12: 1}
    assert factorint(30, limit=2) == {2: 1, 15: 1}
    assert factorint(16, limit=2) == {2: 4}
    assert factorint(124, limit=3) == {2: 2, 31: 1}
    assert factorint(4*31**2, limit=3) == {2: 2, 31: 2}
    p1 = nextprime(2**32)
    p2 = nextprime(2**16)
    p3 = nextprime(p2)
    assert factorint(p1*p2*p3) == {p1: 1, p2: 1, p3: 1}
    assert factorint(13*17*19, limit=15) == {13: 1, 17*19: 1}
    assert factorint(1951*15013*15053, limit=2000) == {225990689: 1, 1951: 1}
    assert factorint(primorial(17) + 1, use_pm1=0) == \
        {long(19026377261): 1, 3467: 1, 277: 1, 105229: 1}
    # when prime b is closer than approx sqrt(8*p) to prime p then they are
    # "close" and have a trivial factorization
    a = nextprime(2**2**8)  # 78 digits
    b = nextprime(a + 2**2**4)
    assert 'Fermat' in capture(lambda: factorint(a*b, verbose=1))

    raises(ValueError, lambda: pollard_rho(4))
    raises(ValueError, lambda: pollard_pm1(3))
    raises(ValueError, lambda: pollard_pm1(10, B=2))
    # verbose coverage
    n = nextprime(2**16)*nextprime(2**17)*nextprime(1901)
    assert 'with primes' in capture(lambda: factorint(n, verbose=1))
    capture(lambda: factorint(nextprime(2**16)*1012, verbose=1))

    n = nextprime(2**17)
    capture(lambda: factorint(n**3, verbose=1))  # perfect power termination
    capture(lambda: factorint(2*n, verbose=1))  # factoring complete msg

    # exceed 1st
    n = nextprime(2**17)
    n *= nextprime(n)
    assert '1000' in capture(lambda: factorint(n, limit=1000, verbose=1))
    n *= nextprime(n)
    assert len(factorint(n)) == 3
    assert len(factorint(n, limit=p1)) == 3
    n *= nextprime(2*n)
    # exceed 2nd
    assert '2001' in capture(lambda: factorint(n, limit=2000, verbose=1))
    assert capture(
        lambda: factorint(n, limit=4000, verbose=1)).count('Pollard') == 2
    # non-prime pm1 result
    n = nextprime(8069)
    n *= nextprime(2*n)*nextprime(2*n, 2)
    capture(lambda: factorint(n, verbose=1))  # non-prime pm1 result
    # factor fermat composite
    p1 = nextprime(2**17)
    p2 = nextprime(2*p1)
    assert factorint((p1*p2**2)**3) == {p1: 3, p2: 6}
    # Test for non integer input
    raises(ValueError, lambda: factorint(4.5))

def test_divisors_and_divisor_count():
    assert divisors(-1) == [1]
    assert divisors(0) == []
    assert divisors(1) == [1]
    assert divisors(2) == [1, 2]
    assert divisors(3) == [1, 3]
    assert divisors(17) == [1, 17]
    assert divisors(10) == [1, 2, 5, 10]
    assert divisors(100) == [1, 2, 4, 5, 10, 20, 25, 50, 100]
    assert divisors(101) == [1, 101]

    assert divisor_count(0) == 0
    assert divisor_count(-1) == 1
    assert divisor_count(1) == 1
    assert divisor_count(6) == 4
    assert divisor_count(12) == 6

    assert divisor_count(180, 3) == divisor_count(180//3)
    assert divisor_count(2*3*5, 7) == 0

def test_issue_6981():
    S = set(divisors(4)).union(set(divisors(Integer(2))))
    assert S == set([1,2,4])


def test_totient():
    assert [totient(k) for k in range(1, 12)] == \
        [1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10]
    assert totient(5005) == 2880
    assert totient(5006) == 2502
    assert totient(5009) == 5008
    assert totient(2**100) == 2**99

    m = Symbol("m", integer=True)
    assert totient(m)
    assert totient(m).subs(m, 3**10) == 3**10 - 3**9
    assert summation(totient(m), (m, 1, 11)) == 42

    n = Symbol("n", integer=True, positive=True)
    assert totient(n).is_integer


def test_divisor_sigma():
    assert [divisor_sigma(k) for k in range(1, 12)] == \
        [1, 3, 4, 7, 6, 12, 8, 15, 13, 18, 12]
    assert [divisor_sigma(k, 2) for k in range(1, 12)] == \
        [1, 5, 10, 21, 26, 50, 50, 85, 91, 130, 122]
    assert divisor_sigma(23450) == 50592
    assert divisor_sigma(23450, 0) == 24
    assert divisor_sigma(23450, 1) == 50592
    assert divisor_sigma(23450, 2) == 730747500
    assert divisor_sigma(23450, 3) == 14666785333344

    m = Symbol("m", integer=True)
    k = Symbol("k", integer=True)
    assert divisor_sigma(m)
    assert divisor_sigma(m, k)
    assert divisor_sigma(m).subs(m, 3**10) == 88573
    assert divisor_sigma(m, k).subs([(m, 3**10), (k, 3)]) == 213810021790597
    assert summation(divisor_sigma(m), (m, 1, 11)) == 99


def test_partitions():
    assert [npartitions(k) for k in range(13)] == \
        [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77]
    assert npartitions(100) == 190569292
    assert npartitions(200) == 3972999029388
    assert npartitions(1000) == 24061467864032622473692149727991
    assert npartitions(2000) == 4720819175619413888601432406799959512200344166
    assert npartitions(10000) % 10**10 == 6916435144
    assert npartitions(100000) % 10**10 == 9421098519


def test_residue():
    assert n_order(2, 13) == 12
    assert [n_order(a, 7) for a in range(1, 7)] == \
           [1, 3, 6, 3, 6, 2]
    assert n_order(5, 17) == 16
    assert n_order(17, 11) == n_order(6, 11)
    assert n_order(101, 119) == 6
    assert n_order(11, (10**50 + 151)**2) == 10000000000000000000000000000000000000000000000030100000000000000000000000000000000000000000000022650
    raises(ValueError, lambda: n_order(6, 9))

    assert is_primitive_root(2, 7) is False
    assert is_primitive_root(3, 8) is False
    assert is_primitive_root(11, 14) is False
    assert is_primitive_root(12, 17) == is_primitive_root(29, 17)
    raises(ValueError, lambda: is_primitive_root(3, 6))

    assert [primitive_root(i) for i in range(2, 31)] == [1, 2, 3, 2, 5, 3, \
       None, 2, 3, 2, None, 2, 3, None, None, 3, 5, 2, None, None, 7, 5, \
       None, 2, 7, 2, None, 2, None]

    for p in primerange(3, 100):
        it = _primitive_root_prime_iter(p)
        assert len(list(it)) == totient(totient(p))
    assert primitive_root(97) == 5
    assert primitive_root(97**2) == 5
    assert primitive_root(40487) == 5
    # note that primitive_root(40487) + 40487 = 40492 is a primitive root
    # of 40487**2, but it is not the smallest
    assert primitive_root(40487**2) == 10
    assert primitive_root(82) == 7
    p = 10**50 + 151
    assert primitive_root(p) == 11
    assert primitive_root(2*p) == 11
    assert primitive_root(p**2) == 11
    raises(ValueError, lambda: primitive_root(-3))

    assert is_quad_residue(3, 7) is False
    assert is_quad_residue(10, 13) is True
    assert is_quad_residue(12364, 139) == is_quad_residue(12364 % 139, 139)
    assert is_quad_residue(207, 251) is True
    assert is_quad_residue(0, 1) is True
    assert is_quad_residue(1, 1) is True
    assert is_quad_residue(0, 2) == is_quad_residue(1, 2) is True
    assert is_quad_residue(1, 4) is True
    assert is_quad_residue(2, 27) is False
    assert is_quad_residue(13122380800, 13604889600) is True
    assert [j for j in range(14) if is_quad_residue(j, 14)] == \
           [0, 1, 2, 4, 7, 8, 9, 11]
    raises(ValueError, lambda: is_quad_residue(1.1, 2))
    raises(ValueError, lambda: is_quad_residue(2, 0))


    assert quadratic_residues(12) == [0, 1, 4, 9]
    assert quadratic_residues(13) == [0, 1, 3, 4, 9, 10, 12]
    assert [len(quadratic_residues(i)) for i in range(1, 20)] == \
      [1, 2, 2, 2, 3, 4, 4, 3, 4, 6, 6, 4, 7, 8, 6, 4, 9, 8, 10]

    assert list(sqrt_mod_iter(6, 2)) == [0]
    assert sqrt_mod(3, 13) == 4
    assert sqrt_mod(3, -13) == 4
    assert sqrt_mod(6, 23) == 11
    assert sqrt_mod(345, 690) == 345

    for p in range(3, 100):
        d = defaultdict(list)
        for i in range(p):
            d[pow(i, 2, p)].append(i)
        for i in range(1, p):
            it = sqrt_mod_iter(i, p)
            v = sqrt_mod(i, p, True)
            if v:
                v = sorted(v)
                assert d[i] == v
            else:
                assert not d[i]

    assert sqrt_mod(9, 27, True) == [3, 6, 12, 15, 21, 24]
    assert sqrt_mod(9, 81, True) == [3, 24, 30, 51, 57, 78]
    assert sqrt_mod(9, 3**5, True) == [3, 78, 84, 159, 165, 240]
    assert sqrt_mod(81, 3**4, True) == [0, 9, 18, 27, 36, 45, 54, 63, 72]
    assert sqrt_mod(81, 3**5, True) == [9, 18, 36, 45, 63, 72, 90, 99, 117,\
            126, 144, 153, 171, 180, 198, 207, 225, 234]
    assert sqrt_mod(81, 3**6, True) == [9, 72, 90, 153, 171, 234, 252, 315,\
            333, 396, 414, 477, 495, 558, 576, 639, 657, 720]
    assert sqrt_mod(81, 3**7, True) == [9, 234, 252, 477, 495, 720, 738, 963,\
            981, 1206, 1224, 1449, 1467, 1692, 1710, 1935, 1953, 2178]

    for a, p in [(26214400, 32768000000), (26214400, 16384000000),
        (262144, 1048576), (87169610025, 163443018796875),
        (22315420166400, 167365651248000000)]:
        assert pow(sqrt_mod(a, p), 2, p) == a

    n = 70
    a, p = 5**2*3**n*2**n, 5**6*3**(n+1)*2**(n+2)
    it = sqrt_mod_iter(a, p)
    for i in range(10):
        assert pow(next(it), 2, p) == a
    a, p = 5**2*3**n*2**n, 5**6*3**(n+1)*2**(n+3)
    it = sqrt_mod_iter(a, p)
    for i in range(2):
        assert pow(next(it), 2, p) == a
    n = 100
    a, p = 5**2*3**n*2**n, 5**6*3**(n+1)*2**(n+1)
    it = sqrt_mod_iter(a, p)
    for i in range(2):
        assert pow(next(it), 2, p) == a

    assert type(next(sqrt_mod_iter(9, 27))) is int
    assert type(next(sqrt_mod_iter(9, 27, ZZ))) is type(ZZ(1))
    assert type(next(sqrt_mod_iter(1, 7, ZZ))) is type(ZZ(1))

    assert is_nthpow_residue(2, 1, 5)
    assert not is_nthpow_residue(2, 2, 5)
    assert is_nthpow_residue(8547, 12, 10007)
    assert nthroot_mod(1801, 11, 2663) == 44
    for a, q, p in [(51922, 2, 203017), (43, 3, 109), (1801, 11, 2663),
          (26118163, 1303, 33333347), (1499, 7, 2663), (595, 6, 2663),
          (1714, 12, 2663), (28477, 9, 33343)]:
        r = nthroot_mod(a, q, p)
        assert pow(r, q, p) == a
    assert nthroot_mod(11, 3, 109) is None

    for p in primerange(5, 100):
        qv = range(3, p, 4)
        for q in qv:
            d = defaultdict(list)
            for i in range(p):
                d[pow(i, q, p)].append(i)
            for a in range(1, p - 1):
                res = nthroot_mod(a, q, p, True)
                if d[a]:
                    assert d[a] == res
                else:
                    assert res is None

    assert legendre_symbol(5, 11) == 1
    assert legendre_symbol(25, 41) == 1
    assert legendre_symbol(67, 101) == -1
    assert legendre_symbol(0, 13) == 0
    assert legendre_symbol(9, 3) == 0
    raises(ValueError, lambda: legendre_symbol(2, 4))

    assert jacobi_symbol(25, 41) == 1
    assert jacobi_symbol(-23, 83) == -1
    assert jacobi_symbol(3, 9) == 0
    assert jacobi_symbol(42, 97) == -1
    assert jacobi_symbol(3, 5) == -1
    assert jacobi_symbol(7, 9) == 1
    assert jacobi_symbol(0, 3) == 0
    assert jacobi_symbol(0, 1) == 1
    assert jacobi_symbol(2, 1) == 1
    assert jacobi_symbol(1, 3) == 1
    raises(ValueError, lambda: jacobi_symbol(3, 8))

    assert mobius(13*7) == 1
    assert mobius(1) == 1
    assert mobius(13*7*5) == -1
    assert mobius(13**2) == 0
    raises(ValueError, lambda: mobius(-3))

    p = Symbol('p', integer=True, positive=True, prime=True)
    x = Symbol('x', positive=True)
    i = Symbol('i', integer=True)
    assert mobius(p) == -1
    raises(TypeError, lambda: mobius(x))
    raises(ValueError, lambda: mobius(i))


def test_hex_pi_nth_digits():
    assert pi_hex_digits(0) == '3243f6a8885a30'
    assert pi_hex_digits(1) == '243f6a8885a308'
    assert pi_hex_digits(10000) == '68ac8fcfb8016c'


def test_crt():
    def mcrt(m, v, r, symmetric=False):
        assert crt(m, v, symmetric)[0] == r
        mm, e, s = crt1(m)
        assert crt2(m, v, mm, e, s, symmetric) == (r, mm)

    mcrt([2, 3, 5], [0, 0, 0], 0)
    mcrt([2, 3, 5], [1, 1, 1], 1)

    mcrt([2, 3, 5], [-1, -1, -1], -1, True)
    mcrt([2, 3, 5], [-1, -1, -1], 2*3*5 - 1, False)

    assert crt([656, 350], [811, 133], symmetric=True) == (-56917, 114800)


def test_binomial_coefficients_list():
    assert binomial_coefficients_list(0) == [1]
    assert binomial_coefficients_list(1) == [1, 1]
    assert binomial_coefficients_list(2) == [1, 2, 1]
    assert binomial_coefficients_list(3) == [1, 3, 3, 1]
    assert binomial_coefficients_list(4) == [1, 4, 6, 4, 1]
    assert binomial_coefficients_list(5) == [1, 5, 10, 10, 5, 1]
    assert binomial_coefficients_list(6) == [1, 6, 15, 20, 15, 6, 1]


def test_binomial_coefficients():
    for n in range(15):
        c = binomial_coefficients(n)
        l = [c[k] for k in sorted(c)]
        assert l == binomial_coefficients_list(n)


def test_multinomial_coefficients():
    assert multinomial_coefficients(1, 1) == {(1,): 1}
    assert multinomial_coefficients(1, 2) == {(2,): 1}
    assert multinomial_coefficients(1, 3) == {(3,): 1}
    assert multinomial_coefficients(2, 0) == {(0, 0): 1}
    assert multinomial_coefficients(2, 1) == {(0, 1): 1, (1, 0): 1}
    assert multinomial_coefficients(2, 2) == {(2, 0): 1, (0, 2): 1, (1, 1): 2}
    assert multinomial_coefficients(2, 3) == {(3, 0): 1, (1, 2): 3, (0, 3): 1,
            (2, 1): 3}
    assert multinomial_coefficients(3, 1) == {(1, 0, 0): 1, (0, 1, 0): 1,
            (0, 0, 1): 1}
    assert multinomial_coefficients(3, 2) == {(0, 1, 1): 2, (0, 0, 2): 1,
            (1, 1, 0): 2, (0, 2, 0): 1, (1, 0, 1): 2, (2, 0, 0): 1}
    mc = multinomial_coefficients(3, 3)
    assert mc == {(2, 1, 0): 3, (0, 3, 0): 1,
            (1, 0, 2): 3, (0, 2, 1): 3, (0, 1, 2): 3, (3, 0, 0): 1,
            (2, 0, 1): 3, (1, 2, 0): 3, (1, 1, 1): 6, (0, 0, 3): 1}
    assert dict(multinomial_coefficients_iterator(2, 0)) == {(0, 0): 1}
    assert dict(
        multinomial_coefficients_iterator(2, 1)) == {(0, 1): 1, (1, 0): 1}
    assert dict(multinomial_coefficients_iterator(2, 2)) == \
        {(2, 0): 1, (0, 2): 1, (1, 1): 2}
    assert dict(multinomial_coefficients_iterator(3, 3)) == mc
    it = multinomial_coefficients_iterator(7, 2)
    assert [next(it) for i in range(4)] == \
        [((2, 0, 0, 0, 0, 0, 0), 1), ((1, 1, 0, 0, 0, 0, 0), 2),
      ((0, 2, 0, 0, 0, 0, 0), 1), ((1, 0, 1, 0, 0, 0, 0), 2)]


def test_issue_4356():
    assert factorint(1030903) == {53: 2, 367: 1}


def test_divisors():
    assert divisors(28) == [1, 2, 4, 7, 14, 28]
    assert [x for x in divisors(3*5*7, 1)] == [1, 3, 5, 15, 7, 21, 35, 105]
    assert divisors(0) == []


def test_divisor_count():
    assert divisor_count(0) == 0
    assert divisor_count(6) == 4


def test_antidivisors():
    assert antidivisors(-1) == []
    assert antidivisors(-3) == [2]
    assert antidivisors(14) == [3, 4, 9]
    assert antidivisors(237) == [2, 5, 6, 11, 19, 25, 43, 95, 158]
    assert antidivisors(12345) == [2, 6, 7, 10, 30, 1646, 3527, 4938, 8230]
    assert antidivisors(393216) == [262144]
    assert sorted(x for x in antidivisors(3*5*7, 1)) == \
        [2, 6, 10, 11, 14, 19, 30, 42, 70]
    assert antidivisors(1) == []


def test_antidivisor_count():
    assert antidivisor_count(0) == 0
    assert antidivisor_count(-1) == 0
    assert antidivisor_count(-4) == 1
    assert antidivisor_count(20) == 3
    assert antidivisor_count(25) == 5
    assert antidivisor_count(38) == 7
    assert antidivisor_count(180) == 6
    assert antidivisor_count(2*3*5) == 3


def test_primorial():
    assert primorial(1) == 2
    assert primorial(1, nth=0) == 1
    assert primorial(2) == 6
    assert primorial(2, nth=0) == 2
    assert primorial(4, nth=0) == 6


def test_smoothness_and_smoothness_p():
    assert smoothness(1) == (1, 1)
    assert smoothness(2**4*3**2) == (3, 16)

    assert smoothness_p(10431, m=1) == \
        (1, [(3, (2, 2, 4)), (19, (1, 5, 5)), (61, (1, 31, 31))])
    assert smoothness_p(10431) == \
        (-1, [(3, (2, 2, 2)), (19, (1, 3, 9)), (61, (1, 5, 5))])
    assert smoothness_p(10431, power=1) == \
        (-1, [(3, (2, 2, 2)), (61, (1, 5, 5)), (19, (1, 3, 9))])
    assert smoothness_p(21477639576571, visual=1) == \
        'p**i=4410317**1 has p-1 B=1787, B-pow=1787\n' + \
        'p**i=4869863**1 has p-1 B=2434931, B-pow=2434931'


def test_visual_factorint():
    assert factorint(1, visual=1) == 1
    forty2 = factorint(42, visual=True)
    assert type(forty2) == Mul
    assert str(forty2) == '2**1*3**1*7**1'
    assert factorint(1, visual=True) is S.One
    no = dict(evaluate=False)
    assert factorint(42**2, visual=True) == Mul(Pow(2, 2, **no),
                                                Pow(3, 2, **no),
                                                Pow(7, 2, **no), **no)
    assert -1 in factorint(-42, visual=True).args


def test_visual_io():
    sm = smoothness_p
    fi = factorint
    # with smoothness_p
    n = 124
    d = fi(n)
    m = fi(d, visual=True)
    t = sm(n)
    s = sm(t)
    for th in [d, s, t, n, m]:
        assert sm(th, visual=True) == s
        assert sm(th, visual=1) == s
    for th in [d, s, t, n, m]:
        assert sm(th, visual=False) == t
    assert [sm(th, visual=None) for th in [d, s, t, n, m]] == [s, d, s, t, t]
    assert [sm(th, visual=2) for th in [d, s, t, n, m]] == [s, d, s, t, t]

    # with factorint
    for th in [d, m, n]:
        assert fi(th, visual=True) == m
        assert fi(th, visual=1) == m
    for th in [d, m, n]:
        assert fi(th, visual=False) == d
    assert [fi(th, visual=None) for th in [d, m, n]] == [m, d, d]
    assert [fi(th, visual=0) for th in [d, m, n]] == [m, d, d]

    # test reevaluation
    no = dict(evaluate=False)
    assert sm({4: 2}, visual=False) == sm(16)
    assert sm(Mul(*[Pow(k, v, **no) for k, v in {4: 2, 2: 6}.items()], **no),
              visual=False) == sm(2**10)

    assert fi({4: 2}, visual=False) == fi(16)
    assert fi(Mul(*[Pow(k, v, **no) for k, v in {4: 2, 2: 6}.items()], **no),
              visual=False) == fi(2**10)


def test_modular():
    assert solve_congruence(*list(zip([3, 4, 2], [12, 35, 17]))) == (1719, 7140)
    assert solve_congruence(*list(zip([3, 4, 2], [12, 6, 17]))) is None
    assert solve_congruence(*list(zip([3, 4, 2], [13, 7, 17]))) == (172, 1547)
    assert solve_congruence(*list(zip([-10, -3, -15], [13, 7, 17]))) == (172, 1547)
    assert solve_congruence(*list(zip([-10, -3, 1, -15], [13, 7, 7, 17]))) is None
    assert solve_congruence(
        *list(zip([-10, -5, 2, -15], [13, 7, 7, 17]))) == (835, 1547)
    assert solve_congruence(
        *list(zip([-10, -5, 2, -15], [13, 7, 14, 17]))) == (2382, 3094)
    assert solve_congruence(
        *list(zip([-10, 2, 2, -15], [13, 7, 14, 17]))) == (2382, 3094)
    assert solve_congruence(*list(zip((1, 1, 2), (3, 2, 4)))) is None
    raises(
        ValueError, lambda: solve_congruence(*list(zip([3, 4, 2], [12.1, 35, 17]))))


def test_search():
    assert 2 in sieve
    assert 2.1 not in sieve
    assert 1 not in sieve
    assert 2**1000 not in sieve
    raises(ValueError, lambda: sieve.search(1))


def test_sieve_slice():
    assert sieve[5] == 11
    assert list(sieve[5:10]) == [sieve[x] for x in range(5, 10)]
    assert list(sieve[5:10:2]) == [sieve[x] for x in range(5, 10, 2)]


def test_continued_fraction():
    raises(ValueError, lambda: cf_p(1, 0, 0))
    raises(ValueError, lambda: cf_p(1, 1, -1))
    assert cf_p(4, 3, 0) == [1, 3]
    assert cf_p(0, 3, 5) == [0, 1, [2, 1, 12, 1, 2, 2]]
    assert cf_p(1, 1, 0) == [1]
    assert cf_p(3, 4, 0) == [0, 1, 3]
    assert cf_p(4, 5, 0) == [0, 1, 4]
    assert cf_p(5, 6, 0) == [0, 1, 5]
    assert cf_p(11, 13, 0) == [0, 1, 5, 2]
    assert cf_p(16, 19, 0) == [0, 1, 5, 3]
    assert cf_p(27, 32, 0) == [0, 1, 5, 2, 2]
    assert cf_p(1, 2, 5) == [[1]]
    assert cf_p(0, 1, 2) == [1, [2]]
    assert cf_p(6, 7, 49) == [1, 1, 6]
    assert cf_p(3796, 1387, 0) == [2, 1, 2, 1, 4]
    assert cf_p(3245, 10000) == [0, 3, 12, 4, 13]
    assert cf_p(1932, 2568) == [0, 1, 3, 26, 2]
    assert cf_p(6589, 2569) == [2, 1, 1, 3, 2, 1, 3, 1, 23]

    def take(iterator, n=7):
        res = []
        for i, t in enumerate(cf_i(iterator)):
            if i >= n:
                break
            res.append(t)
        return res

    assert take(phi) == [1, 1, 1, 1, 1, 1, 1]
    assert take(pi) == [3, 7, 15, 1, 292, 1, 1]

    assert list(cf_i(S(17)/12)) == [1, 2, 2, 2]
    assert list(cf_i(S(-17)/12)) == [-2, 1, 1, 2, 2]

    assert list(cf_c([1, 6, 1, 8])) == [S(1), S(7)/6, S(8)/7, S(71)/62]
    assert list(cf_c([2])) == [S(2)]
    assert list(cf_c([1, 1, 1, 1, 1, 1, 1])) == [S.One, S(2), S(3)/2, S(5)/3,
                                                 S(8)/5, S(13)/8, S(21)/13]
    assert list(cf_c([1, 6, S(-1)/2, 4])) == [S.One, S(7)/6, S(5)/4, S(3)/2]

    assert cf_r([1, 6, 1, 8]) == S(71)/62
    assert cf_r([3]) == S(3)
    assert cf_r([-1, 5, 1, 4]) == S(-24)/29
    assert (cf_r([0, 1, 1, 7, [24, 8]]) - (sqrt(3) + 2)/7).expand() == 0
    assert cf_r([1, 5, 9]) == S(55)/46
    assert (cf_r([[1]]) - (sqrt(5) + 1)/2).expand() == 0


def test_egyptian_fraction():
    def test_equality(r, alg="Greedy"):
        return r == Add(*[Rational(1, i) for i in egyptian_fraction(r, alg)])

    r = random_complex_number(a=0, c=1, b=0, d=0, rational=True)
    assert test_equality(r)

    assert egyptian_fraction(Rational(4, 17)) == [5, 29, 1233, 3039345]
    assert egyptian_fraction(Rational(7, 13), "Greedy") == [2, 26]
    assert egyptian_fraction(Rational(23, 101), "Greedy") == \
        [5, 37, 1438, 2985448, 40108045937720]
    assert egyptian_fraction(Rational(18, 23), "Takenouchi") == \
        [2, 6, 12, 35, 276, 2415]
    assert egyptian_fraction(Rational(5, 6), "Graham Jewett") == \
        [6, 7, 8, 9, 10, 42, 43, 44, 45, 56, 57, 58, 72, 73, 90, 1806, 1807,
         1808, 1892, 1893, 1980, 3192, 3193, 3306, 5256, 3263442, 3263443,
         3267056, 3581556, 10192056, 10650056950806]
    assert egyptian_fraction(Rational(5, 6), "Golomb") == [2, 6, 12, 20, 30]
    assert egyptian_fraction(Rational(5, 121), "Golomb") == [25, 1225, 3577, 7081, 11737]
    raises(ValueError, lambda: egyptian_fraction(Rational(-4, 9)))
    assert egyptian_fraction(Rational(8, 3), "Golomb") == [1, 2, 3, 4, 5, 6, 7,
                                                           14, 574, 2788, 6460,
                                                           11590, 33062, 113820]
    assert egyptian_fraction(Rational(355, 113)) == [1, 2, 3, 4, 5, 6, 7, 8, 9,
                                                     10, 11, 12, 27, 744, 893588,
                                                     1251493536607,
                                                     20361068938197002344405230]


def test_core():
    assert core(35**13, 10) == 42875
    assert core(210**2) == 1
    assert core(7776, 3) == 36
    assert core(10**27, 22) == 10**5
    assert core(537824) == 14
    assert core(1, 6) == 1
