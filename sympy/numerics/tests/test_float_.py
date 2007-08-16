import sys
sys.path.append(".")
import py
from sympy import *
from sympy.numerics import *

def test_bitcount():
    assert bitcount(0) == 0
    assert bitcount(1) == 1
    assert bitcount(7) == 3
    assert bitcount(8) == 4
    assert bitcount(2**100) == 101
    assert bitcount(2**100-1) == 100
    assert bitcount(-(2**100)) == 101
    assert bitcount(-(2**100-1)) == 100

def test_trailing_zeros():
    assert trailing_zeros(0) == 0
    assert trailing_zeros(1) == 0
    assert trailing_zeros(2) == 1
    assert trailing_zeros(7) == 0
    assert trailing_zeros(8) == 3
    assert trailing_zeros(2**100) == 100
    assert trailing_zeros(2**100-1) == 0

def test_round_down():
    assert normalize(0, -4, 4, ROUND_DOWN) == (0, 0)
    assert normalize(0xf0, -4, 4, ROUND_DOWN) == (15, 0)
    assert normalize(0xf1, -4, 4, ROUND_DOWN) == (15, 0)
    assert normalize(0xff, -4, 4, ROUND_DOWN) == (15, 0)
    assert normalize(-0xf0, -4, 4, ROUND_DOWN) == (-15, 0)
    assert normalize(-0xf1, -4, 4, ROUND_DOWN) == (-15, 0)
    assert normalize(-0xff, -4, 4, ROUND_DOWN) == (-15, 0)

def test_round_up():
    assert normalize(0, -4, 4, ROUND_UP) == (0, 0)
    assert normalize(0xf0, -4, 4, ROUND_UP) == (15, 0)
    assert normalize(0xf1, -4, 4, ROUND_UP) == (1, 4)
    assert normalize(0xff, -4, 4, ROUND_UP) == (1, 4)
    assert normalize(-0xf0, -4, 4, ROUND_UP) == (-15, 0)
    assert normalize(-0xf1, -4, 4, ROUND_UP) == (-1, 4)
    assert normalize(-0xff, -4, 4, ROUND_UP) == (-1, 4)

def test_round_floor():
    assert normalize(0, -4, 4, ROUND_FLOOR) == (0, 0)
    assert normalize(0xf0, -4, 4, ROUND_FLOOR) == (15, 0)
    assert normalize(0xf1, -4, 4, ROUND_FLOOR) == (15, 0)
    assert normalize(0xff, -4, 4, ROUND_FLOOR) == (15, 0)
    assert normalize(-0xf0, -4, 4, ROUND_FLOOR) == (-15, 0)
    assert normalize(-0xf1, -4, 4, ROUND_FLOOR) == (-1, 4)
    assert normalize(-0xff, -4, 4, ROUND_FLOOR) == (-1, 4)

def test_round_ceiling():
    assert normalize(0, -4, 4, ROUND_CEILING) == (0, 0)
    assert normalize(0xf0, -4, 4, ROUND_CEILING) == (15, 0)
    assert normalize(0xf1, -4, 4, ROUND_CEILING) == (1, 4)
    assert normalize(0xff, -4, 4, ROUND_CEILING) == (1, 4)
    assert normalize(-0xf0, -4, 4, ROUND_CEILING) == (-15, 0)
    assert normalize(-0xf1, -4, 4, ROUND_CEILING) == (-15, 0)
    assert normalize(-0xff, -4, 4, ROUND_CEILING) == (-15, 0)

def test_round_half_up():
    assert normalize(0, -4, 4, ROUND_HALF_UP) == (0, 0)
    assert normalize(0xf0, -4, 4, ROUND_HALF_UP) == (15, 0)
    assert normalize(0xf7, -4, 4, ROUND_HALF_UP) == (15, 0)
    assert normalize(0xf8, -4, 4, ROUND_HALF_UP) == (1, 4)
    assert normalize(0xf9, -4, 4, ROUND_HALF_UP) == (1, 4)
    assert normalize(0xff, -4, 4, ROUND_HALF_UP) == (1, 4)
    assert normalize(-0xf0, -4, 4, ROUND_HALF_UP) == (-15, 0)
    assert normalize(-0xf7, -4, 4, ROUND_HALF_UP) == (-15, 0)
    assert normalize(-0xf8, -4, 4, ROUND_HALF_UP) == (-1, 4)
    assert normalize(-0xf9, -4, 4, ROUND_HALF_UP) == (-1, 4)
    assert normalize(-0xff, -4, 4, ROUND_HALF_UP) == (-1, 4)

def test_round_half_down():
    assert normalize(0, -4, 4, ROUND_HALF_DOWN) == (0, 0)
    assert normalize(0xf0, -4, 4, ROUND_HALF_DOWN) == (15, 0)
    assert normalize(0xf7, -4, 4, ROUND_HALF_DOWN) == (15, 0)
    assert normalize(0xf8, -4, 4, ROUND_HALF_DOWN) == (15, 0)
    assert normalize(0xf9, -4, 4, ROUND_HALF_DOWN) == (1, 4)
    assert normalize(0xff, -4, 4, ROUND_HALF_DOWN) == (1, 4)
    assert normalize(-0xf0, -4, 4, ROUND_HALF_DOWN) == (-15, 0)
    assert normalize(-0xf7, -4, 4, ROUND_HALF_DOWN) == (-15, 0)
    assert normalize(-0xf8, -4, 4, ROUND_HALF_DOWN) == (-15, 0)
    assert normalize(-0xf9, -4, 4, ROUND_HALF_DOWN) == (-1, 4)
    assert normalize(-0xff, -4, 4, ROUND_HALF_DOWN) == (-1, 4)

def test_round_half_even():
    assert normalize(0, -4, 4, ROUND_HALF_EVEN) == (0, 0)
    assert normalize(0xf0, -4, 4, ROUND_HALF_EVEN) == (15, 0)
    assert normalize(0xf7, -4, 4, ROUND_HALF_EVEN) == (15, 0)
    assert normalize(0xf8, -4, 4, ROUND_HALF_EVEN) == (1, 4)    # 1111.1000 -> 10000.0
    assert normalize(0xf9, -4, 4, ROUND_HALF_EVEN) == (1, 4)    # 1111.1001 -> 10000.0
    assert normalize(0xe8, -4, 4, ROUND_HALF_EVEN) == (7, 1)    # 1110.1000 -> 1110.0
    assert normalize(0xe9, -4, 4, ROUND_HALF_EVEN) == (15, 0)     # 1110.1001 -> 1111.0
    assert normalize(-0xf0, -4, 4, ROUND_HALF_EVEN) == (-15, 0)
    assert normalize(-0xf7, -4, 4, ROUND_HALF_EVEN) == (-15, 0)
    assert normalize(-0xf8, -4, 4, ROUND_HALF_EVEN) == (-1, 4)
    assert normalize(-0xf9, -4, 4, ROUND_HALF_EVEN) == (-1, 4)
    assert normalize(-0xe8, -4, 4, ROUND_HALF_EVEN) == (-7, 1)
    assert normalize(-0xe9, -4, 4, ROUND_HALF_EVEN) == (-15, 0)

def test_cmp():
    import random
    random.seed(123)
    for i in range(100):
        x = (random.random()*1000 - 1000) * 2.0**random.randint(-20, 20)
        y = (random.random()*1000 - 1000) * 2.0**random.randint(-20, 20)
        assert (Float(x) < Float(y)) == (x < y)
        assert (Float(x) > Float(y)) == (x > y)
        assert (Float(x) == Float(y)) == (x == y)
        assert (Float(x) <= Float(y)) == (x <= y)
        assert (Float(x) >= Float(y)) == (x >= y)
        assert Float(x) == Float(x)
    assert Float(4.4408920985006262E-16) < Float(1.7763568394002505E-15)
    assert Float(-4.4408920985006262E-16) > Float(-1.7763568394002505E-15)

def test_almost_equal():
    assert Float(1.2).ae(Float(1.20000001), 1e-7)
    assert not Float(1.2).ae(Float(1.20000001), 1e-9)
    assert not Float(-0.7818314824680298).ae(Float(-0.774695868667929))

def test_int():
    for i in range(-100, 100):
        assert int(Float(i)) == i
    assert int(Float(2**500 + 23)) == 2**500

def test_add():
    assert Float(4) + Float(-70) == -66
    assert Float(1) + Float(1.1)/80 == 1 + 1.1/80
    assert Float((1, 10000000000)) + Float(3) == Float((1, 10000000000))
    assert Float(3) + Float((1, 10000000000)) == Float((1, 10000000000))
    assert Float((1, -10000000000)) + Float(3) == Float(3)
    assert Float(3) + Float((1, -10000000000)) == Float(3)
    assert Float(1) + 1e-15 != 1
    assert Float(1) + 1e-20 == 1

def test_contexts():
    Float.store()
    Float.setprec(100)
    Float.store()
    Float.setprec(200)
    assert Float.getprec() == 200
    Float.revert()
    assert Float.getprec() == 100
    Float.revert()
    assert Float.getprec() == 53
    Float.revert()
    assert Float.getprec() == 53

def test_complex():
    # many more tests needed
    assert 1 + ComplexFloat(2) == 3
    assert not ComplexFloat(2).ae(2+1e-13)
    assert ComplexFloat(2+1e-15j).ae(2)
