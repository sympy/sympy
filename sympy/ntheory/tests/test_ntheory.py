import sys
sys.path.append(".")
import py
from sympy import *
from sympy.ntheory import *

def test_isprime():
    s = Sieve()
    s.extend(100000)
    ps = set(s.primerange(2, 100001))
    for n in range(100001):
        assert (n in ps) == isprime(n)
    assert isprime(179424673)
    # Some Mersenne primes
    assert isprime(2**61 - 1)
    assert isprime(2**89 - 1)
    assert isprime(2**607 - 1)
    assert not isprime(2**601 - 1)

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
    assert nextprime(90) == 97
    assert nextprime(10**40) == (10**40 + 121)
    assert prevprime(3) == 2
    assert prevprime(7) == 5
    assert prevprime(97) == 89
    assert prevprime(10**40) == (10**40 - 17)
    assert list(primerange(2, 7)) == [2, 3, 5]
    assert list(primerange(2, 10)) == [2, 3, 5, 7]
    assert list(primerange(1050, 1100)) == [1051, 1061, \
        1063, 1069, 1087, 1091, 1093, 1097]
    s = Sieve()
    for i in range(30, 2350, 376):
        for j in range(2, 5096, 1139):
            A = list(s.primerange(i, i+j))
            B = list(primerange(i, i+j))
            assert A == B

def test_randprime():
    import random
    random.seed(1234)
    assert randprime(2, 3) == 2
    assert randprime(1, 3) == 2
    assert randprime(3, 5) == 3
    assert py.test.raises(ValueError, 'randprime(20, 22)')
    for a in [100, 300, 500, 250000]:
        for b in [100, 300, 500, 250000]:
            p = randprime(a, a+b)
            assert a <= p < (a+b) and isprime(p)

def test_factor():
    assert trial(1) == []
    assert trial(2) == [(2,1)]
    assert trial(3) == [(3,1)]
    assert trial(4) == [(2,2)]
    assert trial(5) == [(5,1)]
    assert trial(128) == [(2,7)]
    assert trial(720) == [(2,4), (3,2), (5,1)]
    assert factorint(123456) == [(2, 6), (3, 1), (643, 1)]
    assert primefactors(123456) == [2, 3, 643]
    assert factorint(-16) == [(-1, 1), (2, 4)]
    assert factorint(2**(2**6) + 1) == [(274177, 1), (67280421310721, 1)]
