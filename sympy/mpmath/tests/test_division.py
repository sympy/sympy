from sympy.mpmath.lib import *
from sympy.mpmath import mpf

from random import randint, choice, seed

all_modes = [round_floor, round_ceiling, round_down, round_up, round_nearest]

fb = from_bstr
fi = from_int
ff = from_float


def test_div_1_3():
    a = fi(1)
    b = fi(3)
    c = fi(-1)

    # floor rounds down, ceiling rounds up
    assert fdiv(a, b, 7, round_floor)   == fb('0.01010101')
    assert fdiv(a, b, 7, round_ceiling) == fb('0.01010110')
    assert fdiv(a, b, 7, round_down)    == fb('0.01010101')
    assert fdiv(a, b, 7, round_up)      == fb('0.01010110')
    assert fdiv(a, b, 7, round_nearest) == fb('0.01010101')

    # floor rounds up, ceiling rounds down
    assert fdiv(c, b, 7, round_floor)   == fb('-0.01010110')
    assert fdiv(c, b, 7, round_ceiling) == fb('-0.01010101')
    assert fdiv(c, b, 7, round_down)    == fb('-0.01010101')
    assert fdiv(c, b, 7, round_up)      == fb('-0.01010110')
    assert fdiv(c, b, 7, round_nearest) == fb('-0.01010101')

def test_fdivi_1_3():
    a = 1
    b = fi(3)
    c = -1
    assert fdivi(a, b, 7, round_floor)   == fb('0.01010101')
    assert fdivi(a, b, 7, round_ceiling) == fb('0.01010110')
    assert fdivi(a, b, 7, round_down)    == fb('0.01010101')
    assert fdivi(a, b, 7, round_up)      == fb('0.01010110')
    assert fdivi(a, b, 7, round_nearest) == fb('0.01010101')
    assert fdivi(c, b, 7, round_floor)   == fb('-0.01010110')
    assert fdivi(c, b, 7, round_ceiling) == fb('-0.01010101')
    assert fdivi(c, b, 7, round_down)    == fb('-0.01010101')
    assert fdivi(c, b, 7, round_up)      == fb('-0.01010110')
    assert fdivi(c, b, 7, round_nearest) == fb('-0.01010101')


def test_div_300():

    q = fi(1000000)
    a = fi(300499999)    # a/q is a little less than a half-integer
    b = fi(300500000)    # b/q exactly a half-integer
    c = fi(300500001)    # c/q is a little more than a half-integer

    # Check nearest integer rounding (prec=9 as 2**8 < 300 < 2**9)

    assert fdiv(a, q, 9, round_down) == fi(300)
    assert fdiv(b, q, 9, round_down) == fi(300)
    assert fdiv(c, q, 9, round_down) == fi(300)
    assert fdiv(a, q, 9, round_up) == fi(301)
    assert fdiv(b, q, 9, round_up) == fi(301)
    assert fdiv(c, q, 9, round_up) == fi(301)

    # Nearest even integer is down
    assert fdiv(a, q, 9, round_nearest) == fi(300)
    assert fdiv(b, q, 9, round_nearest) == fi(300)
    assert fdiv(c, q, 9, round_nearest) == fi(301)

    # Nearest even integer is up
    a = fi(301499999)
    b = fi(301500000)
    c = fi(301500001)
    assert fdiv(a, q, 9, round_nearest) == fi(301)
    assert fdiv(b, q, 9, round_nearest) == fi(302)
    assert fdiv(c, q, 9, round_nearest) == fi(302)


def test_tight_integer_division():
    # Test that integer division at tightest possible precision is exact
    N = 100
    seed(1)
    for i in range(N):
        a = choice([1, -1]) * randint(1, 1<<randint(10, 100))
        b = choice([1, -1]) * randint(1, 1<<randint(10, 100))
        p = a * b
        width = bitcount(abs(b)) - trailing(b)
        a = fi(a); b = fi(b); p = fi(p)
        for mode in all_modes:
            assert fdiv(p, a, width, mode) == b


def test_epsilon_rounding():
    # Verify that fdiv uses infinite precision; this result will
    # appear to be exactly 0.101 to a near-sighted algorithm

    a = fb('0.101' + ('0'*200) + '1')
    b = fb('1.10101')
    c = fmul(a, b, 250, round_floor) # exact
    assert fdiv(c, b, bitcount(a[1]), round_floor) == a # exact

    assert fdiv(c, b, 2, round_down) == fb('0.10')
    assert fdiv(c, b, 3, round_down) == fb('0.101')
    assert fdiv(c, b, 2, round_up) == fb('0.11')
    assert fdiv(c, b, 3, round_up) == fb('0.110')
    assert fdiv(c, b, 2, round_floor) == fb('0.10')
    assert fdiv(c, b, 3, round_floor) == fb('0.101')
    assert fdiv(c, b, 2, round_ceiling) == fb('0.11')
    assert fdiv(c, b, 3, round_ceiling) == fb('0.110')

    # The same for negative numbers
    a = fb('-0.101' + ('0'*200) + '1')
    b = fb('1.10101')
    c = fmul(a, b, 250, round_floor)
    assert fdiv(c, b, bitcount(a[1]), round_floor) == a

    assert fdiv(c, b, 2, round_down) == fb('-0.10')
    assert fdiv(c, b, 3, round_up) == fb('-0.110')

    # Floor goes up, ceiling goes down
    assert fdiv(c, b, 2, round_floor) == fb('-0.11')
    assert fdiv(c, b, 3, round_floor) == fb('-0.110')
    assert fdiv(c, b, 2, round_ceiling) == fb('-0.10')
    assert fdiv(c, b, 3, round_ceiling) == fb('-0.101')


def test_mod():
    assert mpf(234) % 1 == 0
    assert mpf(-3) % 256 == 253
    assert mpf(0.25) % 23490.5 == 0.25
    assert mpf(0.25) % -23490.5 == -23490.25
    assert mpf(-0.25) % 23490.5 == 23490.25
    assert mpf(-0.25) % -23490.5 == -0.25
    # Check that these cases are handled efficiently
    assert mpf('1e10000000000') % 1 == 0
    assert mpf('1.23e-1000000000') % 1 == mpf('1.23e-1000000000')
    # test __rmod__
    assert 3 % mpf('1.75') == 1.25
