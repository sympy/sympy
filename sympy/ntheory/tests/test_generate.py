from sympy import Sieve, sieve
from sympy.core.compatibility import range

from sympy.ntheory import isprime, totient, randprime, nextprime, prevprime, \
    primerange, primepi, prime, primorial
from sympy.ntheory.generate import cycle_length
from sympy.ntheory.primetest import mr
from sympy.utilities.pytest import raises


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


def test_primorial():
    assert primorial(1) == 2
    assert primorial(1, nth=0) == 1
    assert primorial(2) == 6
    assert primorial(2, nth=0) == 2
    assert primorial(4, nth=0) == 6


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
