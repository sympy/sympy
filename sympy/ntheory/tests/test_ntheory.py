from sympy import Sieve, binomial_coefficients, binomial_coefficients_list, \
        multinomial_coefficients, raises
from sympy.ntheory import isprime, n_order, is_primitive_root, \
    is_quad_residue, legendre_symbol, npartitions, totient, trial, \
    factorint, primefactors, divisors, randprime, nextprime, prevprime, \
    primerange, primepi, prime, pollard_rho
from sympy.ntheory.bbp_pi import pi_hex_digits

from sympy.ntheory.modular import crt, crt1, crt2

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
    s = Sieve()
    assert s[10] == 29

def test_randprime():
    import random
    random.seed(1234)
    assert randprime(2, 3) == 2
    assert randprime(1, 3) == 2
    assert randprime(3, 5) == 3
    raises(ValueError, 'randprime(20, 22)')
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
    assert factorint(5951757) == [(3, 1), (7, 1), (29, 2), (337, 1)]
    assert factorint(64015937) == [(7993, 1), (8009, 1)]
    assert divisors(1) == [1]
    assert divisors(2) == [1, 2]
    assert divisors(3) == [1, 3]
    assert divisors(10) == [1, 2, 5, 10]
    assert divisors(100) == [1, 2, 4, 5, 10, 20, 25, 50, 100]
    assert divisors(101) == [1, 101]
    assert pollard_rho(2**64+1, max_iters=1, seed=1) == 274177
    assert pollard_rho(19) is None

def test_totient():
    assert [totient(k) for k in range(1, 12)] == \
        [1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10]
    assert totient(5005) == 2880
    assert totient(5006) == 2502
    assert totient(5009) == 5008

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
    assert n_order(2,13)==12
    assert [n_order(a,7) for a in range(1,7)]==\
           [1,3,6,3,6,2]
    assert n_order(5,17)==16
    assert n_order(17,11)==n_order(6,11)
    assert n_order(101,119)==6

    assert is_primitive_root(2,7)==False
    assert is_primitive_root(3,8)==False
    assert is_primitive_root(11,14)==False
    assert is_primitive_root(12,17)==is_primitive_root(29,17)

    assert is_quad_residue(3,7)==False
    assert is_quad_residue(10,13)==True
    assert is_quad_residue(12364,139)==is_quad_residue(132,139)
    assert is_quad_residue(207,251)==True

    assert legendre_symbol(5,11)==1
    assert legendre_symbol(25,41)==1
    assert legendre_symbol(67,101)==-1

def test_hex_pi_nth_digits():
    assert pi_hex_digits(0) == '3243f6a8885a30'
    assert pi_hex_digits(1) == '243f6a8885a308'
    assert pi_hex_digits(10000) == '68ac8fcfb8016c'

def test_crt():
    assert crt([2, 3, 5], [0, 0, 0]) == 0
    assert crt([2, 3, 5], [1, 1, 1]) == 1

    assert crt([2, 3, 5], [-1, -1, -1], True) == -1
    assert crt([2, 3, 5], [-1, -1, -1], False) == 2*3*5 - 1


def test_binomial_coefficients_list():
    assert binomial_coefficients_list(0) == [1]
    assert binomial_coefficients_list(1) == [1,1]
    assert binomial_coefficients_list(2) == [1,2,1]
    assert binomial_coefficients_list(3) == [1,3,3,1]
    assert binomial_coefficients_list(4) == [1,4,6,4,1]
    assert binomial_coefficients_list(5) == [1,5,10,10,5,1]
    assert binomial_coefficients_list(6) == [1,6,15,20,15,6,1]

def test_binomial_coefficients():
    for n in range(15):
        c = binomial_coefficients(n)
        l = [c[k] for k in sorted(c)]
        assert l==binomial_coefficients_list(n)

def test_multinomial_coefficients():
    assert multinomial_coefficients(1, 1) == {(1,): 1}
    assert multinomial_coefficients(1, 2) == {(2,): 1}
    assert multinomial_coefficients(1, 3) == {(3,): 1}
    assert multinomial_coefficients(2, 1) == {(0, 1): 1, (1, 0): 1}
    assert multinomial_coefficients(2, 2) == {(2, 0): 1, (0, 2): 1, (1, 1): 2}
    assert multinomial_coefficients(2, 3) == {(3, 0): 1, (1, 2): 3, (0, 3): 1,
            (2, 1): 3}
    assert multinomial_coefficients(3, 1) == {(1, 0, 0): 1, (0, 1, 0): 1,
            (0, 0, 1): 1}
    assert multinomial_coefficients(3, 2) == {(0, 1, 1): 2, (0, 0, 2): 1,
            (1, 1, 0): 2, (0, 2, 0): 1, (1, 0, 1): 2, (2, 0, 0): 1}
    assert multinomial_coefficients(3, 3) == {(2, 1, 0): 3, (0, 3, 0): 1,
            (1, 0, 2): 3, (0, 2, 1): 3, (0, 1, 2): 3, (3, 0, 0): 1,
            (2, 0, 1): 3, (1, 2, 0): 3, (1, 1, 1): 6, (0, 0, 3): 1}
